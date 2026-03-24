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
    load_ledger, get_events, approved_productions, release_label,
)

import h5py


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Config keys that differ between waveform analyses by design and should not
# be flagged as inconsistent.
_WAVEFORM_SPECIFIC_KEYS = frozenset({
    "approximant", "waveform-arguments-dict", "waveform-arguments",
    "frequency-domain-source-model", "time-domain-source-model",
    "waveform-generator", "waveform-generator-class",
    "label", "outdir", "webdir",
    "accounting", "accounting-user", "accounting_user",
    "request-memory", "request-cpus", "request-disk",
    "scheduler-env", "scheduler-analysis-time", "extra-kwargs",
    "nr-file", "NR-file",
})

# These keys are safety-critical if they unexpectedly differ.
_HIGH_PRIORITY_KEYS = frozenset({
    "trigger-time", "trigger_time",
    "detectors", "detectors-list",
    "duration", "sampling-frequency",
    "minimum-frequency", "maximum-frequency",
    "reference-frequency", "channel-dict", "psd-dict",
})


def _read_label_config(fh, key):
    """Return the config dict stored under a PESummary label group."""
    if key not in fh or "config_file" not in fh[key]:
        return {}
    cfg = fh[key]["config_file"]

    def _decode(x):
        return x.decode() if isinstance(x, bytes) else str(x)

    if "keys" in cfg and "vals" in cfg:
        return {_decode(k): _decode(v) for k, v in zip(cfg["keys"], cfg["vals"])}

    result = {}
    for k in cfg:
        try:
            raw = cfg[k][()]
            result[k] = raw.decode() if isinstance(raw, bytes) else str(raw)
        except Exception:
            result[k] = "<unreadable>"
    return result


def _check_config_consistency(hdf5_path, labels):
    """Compare per-label configs in *hdf5_path*.

    Returns a list of dicts, each describing one differing key::

        {"key": str, "values": {label: value}, "high_priority": bool}

    ``labels`` should be the C00-prefixed waveform labels (no Mixed).
    Returns an empty list when fewer than two configs are available.
    """
    try:
        with h5py.File(hdf5_path, "r") as fh:
            configs = {
                label: _read_label_config(fh, label)
                for label in labels
                if label in fh
            }
    except Exception:
        return []

    if len(configs) < 2:
        return []

    all_keys = set()
    for cfg in configs.values():
        all_keys.update(cfg.keys())
    comparable = all_keys - _WAVEFORM_SPECIFIC_KEYS

    diffs = []
    for key in sorted(comparable):
        values = {lbl: cfg.get(key, "<absent>") for lbl, cfg in configs.items()}
        if len(set(values.values())) > 1:
            diffs.append({
                "key": key,
                "values": values,
                "high_priority": key in _HIGH_PRIORITY_KEYS,
            })
    return diffs


# ---------------------------------------------------------------------------
# Rendering helpers
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
    productions = list(approved_productions(event))
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

    # Config consistency — compare non-Mixed waveform labels only
    waveform_labels = [f"C00:{l}" for l in expected_labels if l != "Mixed"]
    info["config_diffs"] = _check_config_consistency(hdf5_path, waveform_labels)

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


def collect_spintaylor_event_info(event_name, event, spintaylor_dir):
    """Collect info for a Stage 1 SpinTaylor-only release file."""
    productions = list(approved_bilby_productions(event))
    if len(productions) != 1:
        return None
    _, settings = productions[0]
    approx = settings.get("waveform", {}).get("approximant", "")
    if "IMRPhenomXPHM" not in approx:
        return None

    label = "IMRPhenomXPHM-SpinTaylor"
    hdf5_path = os.path.join(spintaylor_dir, event_name, "C00-Combined.h5")

    info = {
        "event": event_name,
        "hdf5_path": hdf5_path,
        "hdf5_exists": os.path.isfile(hdf5_path),
        "label": f"C00:{label}",
    }

    if not info["hdf5_exists"]:
        return info

    info["hdf5_size"] = os.path.getsize(hdf5_path)
    info["hdf5_hash"] = sha256(hdf5_path)

    key = f"C00:{label}"
    try:
        with h5py.File(hdf5_path, "r") as fh:
            if key in fh:
                try:
                    info["n_samples"] = fh[key]["posterior_samples"].shape[0]
                except (KeyError, AttributeError):
                    info["n_samples"] = None
            else:
                info["label_missing"] = True
    except Exception as exc:
        info["hdf5_error"] = str(exc)

    fits_name = f"C00-{label}.fits.gz"
    fits_path = os.path.join(spintaylor_dir, event_name, fits_name)
    skymap = {"path": fits_path, "exists": os.path.isfile(fits_path)}
    if skymap["exists"] and os.path.getsize(fits_path) > 0:
        skymap["size"] = os.path.getsize(fits_path)
        skymap["hash"] = sha256(fits_path)
    elif skymap["exists"]:
        skymap["exists"] = False
    info["skymap"] = skymap

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
.warn { color: #9a6700; }
.high-priority { color: #cf222e; font-weight: 600; }
"""


def _truncate(s, n=80):
    return s if len(s) <= n else s[:n] + "\u2026"


def _config_diff_html(diffs):
    if not diffs:
        return '<span class="ok">&#10003; consistent</span>'
    items = []
    for d in diffs:
        vals = "; ".join(
            f"<code>{lbl}</code>={_truncate(v)}"
            for lbl, v in sorted(d["values"].items())
        )
        cls = "high-priority" if d["high_priority"] else "warn"
        items.append(f'<li class="{cls}"><code>{d["key"]}</code>: {vals}</li>')
    return f'<ul>{"" .join(items)}</ul>'


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


def _spintaylor_rows_html(spintaylor_info):
    rows = []
    for info in spintaylor_info:
        event = info["event"]
        if not info["hdf5_exists"]:
            rows.append(
                f'<tr><td><strong>{event}</strong></td>'
                f'<td class="missing" colspan="4">C00-Combined.h5 not found</td></tr>'
            )
            continue
        label = info["label"]
        if info.get("label_missing"):
            label_cell = f'<code>{label}</code> <span class="missing">&#10007; missing</span>'
        elif info.get("hdf5_error"):
            label_cell = f'<span class="missing">{info["hdf5_error"]}</span>'
        elif info.get("n_samples") is None:
            label_cell = f'<code>{label}</code> (unreadable)'
        else:
            label_cell = f'<code>{label}</code>: {info["n_samples"]:,} samples'
        s = info.get("skymap", {})
        if s.get("exists"):
            skymap_cell = '<span class="ok">&#10003;</span> ' + os.path.basename(s["path"])
        else:
            skymap_cell = '<span class="missing">&#10007; missing</span>'
        h = info.get("hdf5_hash", "")
        rows.append(f"""<tr>
  <td><strong>{event}</strong></td>
  <td>{label_cell}</td>
  <td>{human_size(info["hdf5_size"])}</td>
  <td><code>{h}</code></td>
  <td>{skymap_cell}</td>
</tr>""")
    return rows


def render_html(events_info, release_dir, spintaylor_info=None):
    rows = []
    for info in events_info:
        event = info["event"]
        if not info["hdf5_exists"]:
            rows.append(
                f'<tr><td><strong>{event}</strong></td>'
                f'<td class="missing" colspan="5">C00-Combined.h5 not found</td></tr>'
            )
            continue
        h = info.get("hdf5_hash", "")
        rows.append(f"""<tr>
  <td><strong>{event}</strong></td>
  <td>{_label_html(info.get("labels", []))}</td>
  <td>{human_size(info["hdf5_size"])}</td>
  <td><code>{h}</code></td>
  <td>{_skymap_html(info.get("skymaps", []))}</td>
  <td>{_config_diff_html(info.get("config_diffs", []))}</td>
</tr>""")

    spintaylor_section = ""
    if spintaylor_info:
        st_rows = _spintaylor_rows_html(spintaylor_info)
        spintaylor_section = f"""
<h2>Stage 1 &#8212; SpinTaylor-only release</h2>
<table>
  <thead><tr>
    <th>Event</th>
    <th>Label / samples</th>
    <th>HDF5 size</th>
    <th>SHA-256</th>
    <th>Skymap</th>
  </tr></thead>
  <tbody>
{chr(10).join(st_rows)}
  </tbody>
</table>
"""

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
<h2>Stage 2B &#8212; Full combined release</h2>
<table>
  <thead><tr>
    <th>Event</th>
    <th>Waveform labels</th>
    <th>HDF5 size</th>
    <th>SHA-256</th>
    <th>Skymaps</th>
    <th>Config consistency</th>
  </tr></thead>
  <tbody>
{"".join(rows)}
  </tbody>
</table>
{spintaylor_section}</body>
</html>
"""


# ---------------------------------------------------------------------------
# Markdown rendering
# ---------------------------------------------------------------------------

def render_markdown(events_info, release_dir, spintaylor_info=None):
    lines = [
        "# Release report",
        "",
        f"Release directory: `{release_dir}`",
        "",
        "## Stage 2B \u2014 Full combined release",
        "",
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

        diffs = info.get("config_diffs", [])
        if not diffs:
            lines.append("- ✅ Config consistent across waveform analyses")
        else:
            lines.append("- Config inconsistencies between waveform analyses:")
            for d in diffs:
                priority = " **[HIGH PRIORITY]**" if d["high_priority"] else ""
                vals = ", ".join(
                    f"`{lbl}`: {_truncate(v)}"
                    for lbl, v in sorted(d["values"].items())
                )
                lines.append(f"  - `{d['key']}`{priority}: {vals}")

        lines.append("")

    # Stage 1 SpinTaylor section
    if spintaylor_info:
        lines += [
            "## Stage 1 \u2014 SpinTaylor-only release",
            "",
            "| Event | Label | Samples | HDF5 size | SHA-256 (truncated) | Skymap |",
            "|-------|-------|---------|-----------|---------------------|--------|",
        ]
        for info in spintaylor_info:
            event = info["event"]
            if not info["hdf5_exists"]:
                lines.append(f"| {event} | \u274c HDF5 missing | \u2014 | \u2014 | \u2014 | \u2014 |")
                continue
            label = info["label"]
            if info.get("label_missing"):
                label_cell = f"`{label}` \u274c"
                count_cell = "\u2014"
            elif info.get("n_samples") is None:
                label_cell = f"`{label}`"
                count_cell = "\u2014"
            else:
                label_cell = f"`{label}`"
                count_cell = f"{info['n_samples']:,}"
            h = info.get("hdf5_hash", "")
            s = info.get("skymap", {})
            skymap_cell = "\u2705" if s.get("exists") else "\u274c"
            lines.append(
                f"| {event} | {label_cell} | {count_cell} | "
                f"{human_size(info['hdf5_size'])} | `{h[:20]}\u2026` | {skymap_cell} |"
            )

        lines += ["", "### Stage 1 file hashes", "", "```"]
        for info in spintaylor_info:
            if not info.get("hdf5_exists") or "hdf5_hash" not in info:
                continue
            lines.append(f"{info['hdf5_hash']}  {info['hdf5_path']}")
            s = info.get("skymap", {})
            if s.get("exists") and "hash" in s:
                lines.append(f"{s['hash']}  {s['path']}")
        lines += ["```", ""]

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate release summary report (HTML + Markdown)")
    parser.add_argument("--ledger", default=".asimov/ledger.yml")
    parser.add_argument("--release-dir", default="releases/staging")
    parser.add_argument("--spintaylor-dir", default=None,
                        help="Stage 1 SpinTaylor release directory "
                             "(default: {release-dir}/spintaylor if it exists)")
    parser.add_argument("--output-dir", default=None,
                        help="Where to write reports (default: --release-dir)")
    args = parser.parse_args()

    output_dir = args.output_dir or args.release_dir
    os.makedirs(output_dir, exist_ok=True)

    spintaylor_dir = args.spintaylor_dir
    if spintaylor_dir is None:
        auto = os.path.join(args.release_dir, "spintaylor")
        if os.path.isdir(auto):
            spintaylor_dir = auto

    data = load_ledger(args.ledger)
    events = get_events(data)

    events_info = []
    for name, event in sorted(events.items()):
        info = collect_event_info(name, event, args.release_dir)
        if info is not None:
            events_info.append(info)

    spintaylor_info = []
    if spintaylor_dir:
        print(f"Including Stage 1 SpinTaylor release from: {spintaylor_dir}")
        for name, event in sorted(events.items()):
            info = collect_spintaylor_event_info(name, event, spintaylor_dir)
            if info is not None:
                spintaylor_info.append(info)

    html_path = os.path.join(output_dir, "release_report.html")
    md_path = os.path.join(output_dir, "release_report.md")

    with open(html_path, "w") as f:
        f.write(render_html(events_info, args.release_dir, spintaylor_info or None))
    print(f"Written: {html_path}")

    with open(md_path, "w") as f:
        f.write(render_markdown(events_info, args.release_dir, spintaylor_info or None))
    print(f"Written: {md_path}")


if __name__ == "__main__":
    main()
