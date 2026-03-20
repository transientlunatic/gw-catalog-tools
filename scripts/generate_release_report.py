#!/usr/bin/env python3
"""Generate HTML and Markdown reports summarising the contents of a release.

For each event with APPROVED bilby productions the report shows:
  - Path to the combined HDF5 file
  - Waveform labels present and posterior sample count per label
  - File size and SHA-256 hash
  - Per-label skymap FITS file, size, and hash

Usage:
    python scripts/generate_release_report.py [--ledger .asimov/ledger.yml]
                                               [--release-dir releases/staging]
                                               [--output-dir releases/staging]
"""

import argparse
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from ledger import (
    load_ledger, get_events, approved_bilby_productions, release_label,
)

import h5py


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def human_size(n):
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"


def _expected_labels(productions):
    labels = [
        release_label(settings["waveform"]["approximant"])
        for _, settings in productions
    ]
    if len(labels) > 1:
        labels.append("Mixed")
    return labels


# ---------------------------------------------------------------------------
# Data collection
# ---------------------------------------------------------------------------

def collect_event_info(event_name, event, release_dir):
    productions = list(approved_bilby_productions(event))
    if not productions:
        return None

    expected_labels = _expected_labels(productions)
    hdf5_path = os.path.join(release_dir, event_name, "C00-Combined.h5")

    info = {
        "event": event_name,
        "hdf5_path": hdf5_path,
        "hdf5_exists": os.path.isfile(hdf5_path),
        "labels": [],
        "skymaps": [],
    }

    if not info["hdf5_exists"]:
        return info

    info["hdf5_size"] = os.path.getsize(hdf5_path)
    info["hdf5_hash"] = sha256(hdf5_path)

    try:
        with h5py.File(hdf5_path, "r") as fh:
            for label in expected_labels:
                key = f"C00:{label}"
                if key in fh:
                    try:
                        n = fh[key]["posterior_samples"].shape[0]
                    except (KeyError, AttributeError):
                        n = None
                    info["labels"].append({"label": key, "n_samples": n})
                else:
                    info["labels"].append({"label": key, "n_samples": None, "missing": True})
    except Exception as exc:
        info["hdf5_error"] = str(exc)

    for label in expected_labels:
        fits_name = f"C00-{label.replace(':', '-')}.fits.gz"
        fits_path = os.path.join(release_dir, event_name, fits_name)
        skymap = {"label": label, "path": fits_path, "exists": os.path.isfile(fits_path)}
        if skymap["exists"] and os.path.getsize(fits_path) > 0:
            skymap["size"] = os.path.getsize(fits_path)
            skymap["hash"] = sha256(fits_path)
        elif skymap["exists"]:
            skymap["exists"] = False  # treat empty as missing
        info["skymaps"].append(skymap)

    return info


# ---------------------------------------------------------------------------
# HTML rendering
# ---------------------------------------------------------------------------

_HTML_STYLE = """
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
       max-width: 1200px; margin: 2em auto; padding: 0 1em; color: #333; }
h1 { border-bottom: 2px solid #e1e4e8; padding-bottom: 0.3em; }
table { border-collapse: collapse; width: 100%; margin-bottom: 2em; }
th, td { border: 1px solid #d0d7de; padding: 0.5em 0.75em; vertical-align: top; }
th { background: #f6f8fa; font-weight: 600; }
tr:nth-child(even) { background: #f6f8fa; }
ul { margin: 0; padding-left: 1.4em; }
code { background: #eef0f2; padding: 0.1em 0.35em; border-radius: 3px;
       font-size: 0.9em; word-break: break-all; }
.missing { color: #cf222e; }
.ok { color: #1a7f37; }
"""


def _label_html(labels):
    items = []
    for l in labels:
        n = l.get("n_samples")
        if l.get("missing"):
            items.append(f'<li><code>{l["label"]}</code> <span class="missing">✗ missing</span></li>')
        elif n is None:
            items.append(f'<li><code>{l["label"]}</code> (unreadable)</li>')
        else:
            items.append(f'<li><code>{l["label"]}</code>: {n:,} samples</li>')
    return f'<ul>{"".join(items)}</ul>'


def _skymap_html(skymaps):
    items = []
    for s in skymaps:
        if s["exists"]:
            items.append(f'<li><span class="ok">✓</span> {s["label"]}</li>')
        else:
            items.append(f'<li><span class="missing">✗</span> {s["label"]}</li>')
    return f'<ul>{"".join(items)}</ul>'


def render_html(events_info, release_dir):
    rows = []
    for info in events_info:
        event = info["event"]
        if not info["hdf5_exists"]:
            rows.append(
                f'<tr><td><strong>{event}</strong></td>'
                f'<td class="missing" colspan="4">C00-Combined.h5 not found</td></tr>'
            )
            continue
        h = info.get("hdf5_hash", "")
        rows.append(f"""<tr>
  <td><strong>{event}</strong></td>
  <td>{_label_html(info.get("labels", []))}</td>
  <td>{human_size(info["hdf5_size"])}</td>
  <td><code>{h}</code></td>
  <td>{_skymap_html(info.get("skymaps", []))}</td>
</tr>""")

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Release report — {release_dir}</title>
  <style>{_HTML_STYLE}</style>
</head>
<body>
<h1>Release report</h1>
<p>Release directory: <code>{release_dir}</code></p>
<table>
  <thead><tr>
    <th>Event</th>
    <th>Waveform labels</th>
    <th>HDF5 size</th>
    <th>SHA-256</th>
    <th>Skymaps</th>
  </tr></thead>
  <tbody>
{"".join(rows)}
  </tbody>
</table>
</body>
</html>
"""


# ---------------------------------------------------------------------------
# Markdown rendering
# ---------------------------------------------------------------------------

def render_markdown(events_info, release_dir):
    lines = [
        f"# Release report",
        f"",
        f"Release directory: `{release_dir}`",
        f"",
    ]

    # Summary table
    lines += [
        "| Event | Labels | Samples | HDF5 size | SHA-256 (truncated) |",
        "|-------|--------|---------|-----------|---------------------|",
    ]
    for info in events_info:
        event = info["event"]
        if not info["hdf5_exists"]:
            lines.append(f"| {event} | ❌ HDF5 missing | — | — | — |")
            continue
        label_strs = []
        for l in info.get("labels", []):
            if l.get("missing"):
                label_strs.append(f"`{l['label']}` ❌")
            else:
                label_strs.append(f"`{l['label']}`")
        labels = "<br>".join(label_strs)
        counts = "<br>".join(
            f"{l['n_samples']:,}" if l.get("n_samples") is not None and not l.get("missing") else "—"
            for l in info.get("labels", [])
        )
        h = info.get("hdf5_hash", "")
        lines.append(
            f"| {event} | {labels} | {counts} | "
            f"{human_size(info['hdf5_size'])} | `{h[:20]}…` |"
        )

    # Hash block
    lines += ["", "## File hashes", "", "```"]
    for info in events_info:
        if not info.get("hdf5_exists") or "hdf5_hash" not in info:
            continue
        lines.append(f"{info['hdf5_hash']}  {info['hdf5_path']}")
        for s in info.get("skymaps", []):
            if s.get("exists") and "hash" in s:
                lines.append(f"{s['hash']}  {s['path']}")
    lines += ["```", ""]

    # Per-event details
    lines += ["## Per-event details", ""]
    for info in events_info:
        event = info["event"]
        lines.append(f"### {event}")
        lines.append("")
        if not info["hdf5_exists"]:
            lines.append(f"- ❌ `{info['hdf5_path']}` not found")
            lines.append("")
            continue

        lines.append(f"**HDF5**: `{info['hdf5_path']}`  ")
        lines.append(f"Size: {human_size(info['hdf5_size'])}  ")
        lines.append(f"SHA-256: `{info.get('hdf5_hash', 'N/A')}`")
        lines.append("")

        for l in info.get("labels", []):
            if l.get("missing"):
                lines.append(f"- ❌ `{l['label']}` not found in HDF5")
            elif l.get("n_samples") is None:
                lines.append(f"- ⚠️ `{l['label']}` unreadable")
            else:
                lines.append(f"- ✅ `{l['label']}`: {l['n_samples']:,} samples")

        lines.append("")
        for s in info.get("skymaps", []):
            if s.get("exists") and "hash" in s:
                lines.append(f"- ✅ Skymap `{s['label']}`: `{s['path']}`")
                lines.append(f"  - Size: {human_size(s['size'])}, SHA-256: `{s['hash']}`")
            else:
                lines.append(f"- ❌ Skymap `{s['label']}`: missing")
        lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate release summary report (HTML + Markdown)")
    parser.add_argument("--ledger", default=".asimov/ledger.yml")
    parser.add_argument("--release-dir", default="releases/staging")
    parser.add_argument("--output-dir", default=None,
                        help="Where to write reports (default: --release-dir)")
    args = parser.parse_args()

    output_dir = args.output_dir or args.release_dir
    os.makedirs(output_dir, exist_ok=True)

    data = load_ledger(args.ledger)
    events = get_events(data)

    events_info = []
    for name, event in sorted(events.items()):
        info = collect_event_info(name, event, args.release_dir)
        if info is not None:
            events_info.append(info)

    html_path = os.path.join(output_dir, "release_report.html")
    md_path = os.path.join(output_dir, "release_report.md")

    with open(html_path, "w") as f:
        f.write(render_html(events_info, args.release_dir))
    print(f"Written: {html_path}")

    with open(md_path, "w") as f:
        f.write(render_markdown(events_info, args.release_dir))
    print(f"Written: {md_path}")


if __name__ == "__main__":
    main()
