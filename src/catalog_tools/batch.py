"""
Batch review: run the checks over a tree of results and summarise them.

Results are laid out as ``root / event / analysis / posterior_samples.h5``.
:func:`review_file` records the checks on one file, :func:`review_tree`
collects those into the JSON written by ``catalog-tools batch``, and
:func:`generate_markdown` renders that JSON as a table for
``catalog-tools report``.
"""

import fnmatch
from datetime import datetime
from pathlib import Path

from .review import open_pesummary_samples, read_run_info, run_all_checks

DEFAULT_EXCLUDE = ("generate-psd",)

STATUS_RANK = {"ERROR": 4, "FAIL": 3, "WARNING": 2, "INFO": 1, "PASS": 0}

STATUS_TEXT = {
    "PASS": "PASS", "INFO": "INFO", "WARNING": "WARN", "FAIL": "FAIL",
    "ERROR": "ERR", "SKIP": "N/A", "N/A": "N/A",
}

STATUS_EMOJI = {
    "PASS": "✅", "INFO": "🔍", "WARNING": "⚠️", "FAIL": "❌",
    "ERROR": "⛔", "SKIP": "N/A", "N/A": "N/A",
}

# The order of the columns of the report
TEST_ORDER = (
    "effective_sample_size",
    "chirp_mass_railing",
    "distance_railing",
    "snr_negative",
    "spin_prior_comparison",
    "multimodality",
    "diptest_unimodality",
    "parameter_consistency",
    "likelihood_sanity",
    "network_snr",
    "snr_plausibility",
    "sky_localization",
)

TEST_LABELS = {
    "effective_sample_size": "Samples",
    "chirp_mass_railing": "Chirp Mass Railing",
    "distance_railing": "Distance Railing",
    "snr_negative": "SNR Negative",
    "spin_prior_comparison": "Spin Prior",
    "multimodality": "Multimodality",
    "diptest_unimodality": "Dip Test",
    "parameter_consistency": "Param Consistency",
    "likelihood_sanity": "Log Likelihood",
    "network_snr": "Network SNR",
    "snr_plausibility": "SNR Plausibility",
    "sky_localization": "Sky Localization",
}


def worst(statuses):
    """The most severe of the statuses; statuses with no rank (SKIP) count as N/A."""
    ranked = [s for s in statuses if s in STATUS_RANK]
    return max(ranked, key=STATUS_RANK.get) if ranked else "N/A"


def discover_files(root, analysis_patterns=None, exclude=DEFAULT_EXCLUDE):
    """
    Find the ``posterior_samples.h5`` files under ``root``.

    Yields ``(event, analysis, path)``.  ``analysis_patterns`` are names or
    globs (``"bilby-*"``) the analysis directory must match, if any are given;
    analysis directories named in ``exclude`` are always skipped.
    """
    root = Path(root)
    if not root.is_dir():
        raise FileNotFoundError(f"Root directory not found: {root}")
    for event_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        for analysis_dir in sorted(p for p in event_dir.iterdir() if p.is_dir()):
            if analysis_dir.name in exclude:
                continue
            if analysis_patterns and not any(fnmatch.fnmatch(analysis_dir.name, p)
                                             for p in analysis_patterns):
                continue
            path = analysis_dir / "posterior_samples.h5"
            if path.exists():
                yield event_dir.name, analysis_dir.name, path


def parse_path_ids(path):
    """``(event, analysis)`` from ``.../event/analysis/posterior_samples.h5``."""
    path = Path(path).resolve()
    return path.parent.parent.name, path.parent.name


def review_file(path):
    """
    Run every check on one file.

    Returns a record with the file, a timestamp, the number of samples, the
    ``load_error`` if the file couldn't be read, ``tests`` (``{name: {status,
    notes}}``) and the ``overall_status``, the worst of the tests.
    """
    record = {
        "file": str(path),
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "n_samples": None,
        "load_error": None,
        "tests": {},
        "overall_status": "N/A",
    }
    try:
        posterior, prior = open_pesummary_samples(path)
        info = read_run_info(path)
    except Exception as error:
        record["load_error"] = str(error)
        return record
    record["n_samples"] = len(posterior)

    try:
        results = run_all_checks(posterior, prior, info)
    except Exception as error:
        record["tests"]["run_all_checks"] = {"status": "ERROR", "notes": str(error)}
    else:
        for name, result in results.items():
            quiet = result["status"] in ("PASS", "INFO", "SKIP")
            record["tests"][name] = {"status": result["status"],
                                     "notes": "" if quiet else result["message"]}
    record["overall_status"] = worst(t["status"] for t in record["tests"].values())
    return record


def review_tree(files, meta=None):
    """
    Review ``files``, an iterable of ``(event, analysis, path)``.

    Returns the data written to JSON: ``{"meta": ..., "events": {event:
    {analysis: record}}}``.
    """
    events = {}
    for event, analysis, path in files:
        events.setdefault(event, {})[analysis] = review_file(path)
    return {"meta": {"generated": datetime.now().strftime("%Y-%m-%d %H:%M:%S"), **(meta or {})},
            "events": events}


def generate_markdown(data, emoji=True, ignore=None):
    """
    A markdown summary of the data from :func:`review_tree`.

    ``ignore`` lists tests to leave out of the table; they remain in the JSON.
    """
    symbols = STATUS_EMOJI if emoji else STATUS_TEXT
    meta = data.get("meta", {})
    lines = ["# Automated PE Review"]
    if meta.get("generated"):
        lines.append(f"\nGenerated: {meta['generated']}")
    if meta.get("root"):
        lines.append(f"Root: `{meta['root']}`")
    if meta.get("analysis_patterns"):
        lines.append("Analysis patterns: " + ", ".join(f"`{p}`" for p in meta["analysis_patterns"]))
    if meta.get("exclude"):
        lines.append("Excluded dirs: " + ", ".join(f"`{e}`" for e in meta["exclude"]))
    lines.append("")

    rows = [(event, analysis, record)
            for event, analyses in sorted(data.get("events", {}).items())
            for analysis, record in sorted(analyses.items())]
    if not rows:
        return "\n".join(lines + ["_No results found._"])

    railing = [(event, analysis, record["tests"][name])
               for event, analysis, record in rows if not record.get("load_error")
               for name in ("chirp_mass_railing", "distance_railing")
               if record["tests"].get(name, {}).get("status") not in (None, "PASS", "SKIP", "N/A")]
    lines.append("## Railing summary\n")
    if railing:
        lines += ["The following analyses have railing results; review the prior bounds:", ""]
        lines += [f"- **{event}** / `{analysis}`" + (f" — {t['notes']}" if t.get("notes") else "")
                  for event, analysis, t in railing]
    else:
        lines.append("_No railing issues detected._")
    lines.append("")

    ignore = set(ignore or [])
    present = {name for _, _, record in rows for name in record.get("tests", {})}
    columns = [t for t in TEST_ORDER if t in present and t not in ignore]
    columns += sorted(present - set(TEST_ORDER) - ignore)
    header = ["Event", "Analysis", "N samples", "Overall"]
    header += [TEST_LABELS.get(t, t.replace("_", " ").title()) for t in columns] + ["Notes"]
    lines.append("| " + " | ".join(header) + " |")
    lines.append("| " + " | ".join(["---"] * len(header)) + " |")

    for event, analysis, record in rows:
        if record.get("load_error"):
            cells = [event, analysis, "—", symbols["ERROR"]] + ["—"] * len(columns)
            cells.append(f"Load error: {record['load_error']}")
        else:
            n = record.get("n_samples")
            overall = record.get("overall_status", "N/A")
            cells = [event, analysis, f"{n:,}" if n is not None else "—", symbols.get(overall, overall)]
            notes = []
            for t in columns:
                result = record["tests"].get(t, {})
                status = result.get("status", "N/A")
                cells.append(symbols.get(status, status))
                if result.get("notes"):
                    notes.append(result["notes"])
            cells.append("; ".join(notes))
        lines.append("| " + " | ".join(cells) + " |")

    counts = {s: sum(1 for _, _, r in rows if r.get("overall_status") == s)
              for s in ("PASS", "INFO", "WARNING", "FAIL")}
    errors = sum(1 for _, _, r in rows if r.get("load_error"))
    lines += ["", f"**Total runs:** {len(rows)} | PASS: {counts['PASS']} | INFO: {counts['INFO']} | "
                  f"WARNING: {counts['WARNING']} | FAIL: {counts['FAIL']} | ERROR: {errors}"]
    return "\n".join(lines)
