#!/usr/bin/env python3
"""
Standalone batch PE review script.

Two modes:

  run     Discover events under a results root, run all quality checks, write JSON.
  report  Read a JSON produced by 'run' and write a markdown summary table.

Examples
--------
# Discover ALL analyses under a root (auto-skips generate-psd):
python test_automated_review_standalone.py run \\
    --root /home/pe.o4/GWTC6/project/results \\
    --output review.json --markdown review.md

# Match specific analyses (one or more names/patterns):
python test_automated_review_standalone.py run \\
    --root /home/pe.o4/GWTC6/project/results \\
    --analysis bilby-IMRPhenomXPHM-SpinTaylor bilby-IMRPhenomXPNR bilby-NRSur7dq4 \\
    --output review.json --markdown review.md

# Glob-match multiple analyses per event:
python test_automated_review_standalone.py run \\
    --root /home/pe.o4/GWTC6/project/results \\
    --analysis "bilby-*" \\
    --output review.json

# Skip extra directories on top of the default (generate-psd):
python test_automated_review_standalone.py run \\
    --root /home/pe.o4/GWTC6/project/results \\
    --exclude generate-psd some-other-dir \\
    --output review.json

# Test explicit H5 files (bypass auto-discovery):
python test_automated_review_standalone.py run \\
    path/to/S250202cu/bilby-IMRPhenomXPHM-SpinTaylor/posterior_samples.h5 \\
    --output review.json

# Generate a markdown report from saved JSON:
python test_automated_review_standalone.py report review.json --output review.md
python test_automated_review_standalone.py report review.json  # stdout
"""

import sys
import fnmatch
import json
import argparse
from datetime import datetime
from pathlib import Path

import numpy as np
import h5py as h5
import peconfigurator.auxiliary.run_analyzer as run_analyzer

from diptest_core import hartigan_diptest

# P-value threshold for unimodality: distributions with pval > this are
# considered unimodal under Hartigan's dip test.
UNIMODAL_PVAL_THRESHOLD = 0.01

# Hyperparameters for the railing check (passed to run_analyzer.railing_check).
RAILING_NBIN = 50
RAILING_TOLERANCE = 2.0


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def open_pesummary_samples(metafile):
    """Open and extract posterior and prior samples from a PESummary HDF5 file."""
    with h5.File(metafile, 'r') as metafile_h:
        analyses = list(metafile_h.keys())
        for skip in ("history", "version"):
            if skip in analyses:
                analyses.remove(skip)

        key = analyses[0]

        posterior = np.array(metafile_h[key]['posterior_samples']) # type: ignore

        prior = {}
        if 'priors' in metafile_h[key] and 'samples' in metafile_h[key]['priors']: # type: ignore
            prior = {
                name: np.array(data)
                for name, data in metafile_h[key]['priors']['samples'].items() # type: ignore
            }

    return posterior, prior


# ---------------------------------------------------------------------------
# Individual test functions
# Each returns a raw dict; _collapse() normalises it to {status, notes}.
# ---------------------------------------------------------------------------

def test_spin_prior_comparison(posterior, prior):
    """Compare spin parameter posteriors to priors."""
    results = {}

    for param in ['a_1', 'a_2', 'chi_eff', 'chi_p']:
        if param not in posterior.dtype.names:
            continue
        if param not in prior or len(prior[param]) == 0:
            continue

        post_data = posterior[param]
        prior_data = prior[param]
        status = "PASS"
        warnings = []

        if param in ['a_1', 'a_2', 'chi_p']:
            frac_high_post = np.sum(post_data > 0.7) / len(post_data)
            frac_high_prior = np.sum(prior_data > 0.7) / len(prior_data)
            if frac_high_post > frac_high_prior + 0.1:
                status = "INTERESTING"
                label = "high spin" if param != 'chi_p' else "high precession"
                warnings.append(
                    f"{param}: posterior {frac_high_post*100:.1f}% > 0.7 "
                    f"({label}), prior {frac_high_prior*100:.1f}%"
                )
            results[param] = {
                'status': status,
                'frac_high_post': frac_high_post,
                'frac_high_prior': frac_high_prior,
                'warnings': warnings,
            }
        else:  # chi_eff
            post_upper = np.sum(post_data > 0.9) / len(post_data)
            post_lower = np.sum(post_data < -0.9) / len(post_data)
            prior_upper = np.sum(prior_data > 0.9) / len(prior_data)
            prior_lower = np.sum(prior_data < -0.9) / len(prior_data)
            if post_upper > prior_upper + 0.1:
                status = "INTERESTING"
                warnings.append(
                    f"chi_eff: posterior {post_upper*100:.1f}% near +1, "
                    f"prior {prior_upper*100:.1f}%"
                )
            if post_lower > prior_lower + 0.1:
                status = "INTERESTING"
                warnings.append(
                    f"chi_eff: posterior {post_lower*100:.1f}% near -1, "
                    f"prior {prior_lower*100:.1f}%"
                )
            results[param] = {
                'status': status,
                'post_upper': post_upper,
                'post_lower': post_lower,
                'prior_upper': prior_upper,
                'prior_lower': prior_lower,
                'warnings': warnings,
            }

    return results


def test_sample_count(posterior):
    """Test that sufficient samples were collected."""
    n_samples = len(posterior)
    fail_threshold = 2000
    warning_threshold = 5000

    if n_samples < fail_threshold:
        status = "FAIL"
    elif n_samples < warning_threshold:
        status = "WARNING"
    else:
        status = "PASS"

    return {
        'n_samples': n_samples,
        'fail_threshold': fail_threshold,
        'warning_threshold': warning_threshold,
        'status': status,
    }


def test_multimodality(posterior):
    """Check for multimodality in key parameters."""
    results = {}

    for param in ['mass_ratio', 'chirp_mass', 'luminosity_distance']:
        if param not in posterior.dtype.names:
            continue

        data = posterior[param]
        sorted_data = np.sort(data)
        diffs = np.diff(sorted_data)
        median_diff = np.median(diffs)
        large_gaps = diffs > 5 * median_diff

        if np.any(large_gaps):
            n_gaps = int(np.sum(large_gaps))
            max_gap = float(np.max(diffs[large_gaps]))
            gap_fraction = float(max_gap / (np.max(data) - np.min(data)))

            if gap_fraction > 0.05:
                status = "WARNING"
            else:
                status = "PASS"

            results[param] = {
                'status': status,
                'n_large_gaps': n_gaps,
                'max_gap_fraction': gap_fraction,
            }
        else:
            results[param] = {'status': 'PASS'}

    return results


def test_parameter_consistency(posterior):
    """Test derived parameter consistency."""
    results = {}

    if all(p in posterior.dtype.names for p in ['total_mass', 'mass_1', 'mass_2']):
        calculated = posterior['mass_1'] + posterior['mass_2']
        actual = posterior['total_mass']
        rel_error = np.abs(actual - calculated) / (np.abs(actual) + 1e-10)
        max_error = float(np.max(rel_error))
        status = "PASS" if max_error < 1e-3 else "FAIL"
        results['total_mass'] = {'status': status, 'max_error': max_error}

    if all(p in posterior.dtype.names for p in ['mass_ratio', 'mass_1', 'mass_2']):
        calculated = posterior['mass_2'] / posterior['mass_1']
        actual = posterior['mass_ratio']
        rel_error = np.abs(actual - calculated) / (np.abs(actual) + 1e-10)
        max_error = float(np.max(rel_error))
        status = "PASS" if max_error < 1e-3 else "FAIL"
        results['mass_ratio'] = {'status': status, 'max_error': max_error}

    if all(p in posterior.dtype.names for p in ['chirp_mass', 'mass_1', 'mass_2']):
        m1 = posterior['mass_1']
        m2 = posterior['mass_2']
        calculated = (m1 * m2) ** (3 / 5) / (m1 + m2) ** (1 / 5)
        actual = posterior['chirp_mass']
        rel_error = np.abs(actual - calculated) / (np.abs(actual) + 1e-10)
        max_error = float(np.max(rel_error))
        status = "PASS" if max_error < 1e-3 else "FAIL"
        results['chirp_mass'] = {'status': status, 'max_error': max_error}

    return results


def test_diptest_unimodality(posterior, unimodal_pval_threshold=UNIMODAL_PVAL_THRESHOLD):
    """Hartigan's dip test for unimodality on key 1-D posterior parameters."""
    results = {}

    for param in ['mass_ratio', 'chirp_mass', 'luminosity_distance']:
        if param not in posterior.dtype.names:
            continue

        data = posterior[param].astype(float)
        if len(data) < 4:
            results[param] = {'status': 'N/A', 'dip': None, 'pval': None}
            continue

        dip, pval = hartigan_diptest(data)
        status = 'PASS' if pval > unimodal_pval_threshold else 'WARNING'
        results[param] = {'status': status, 'dip': float(dip), 'pval': float(pval)}

    return results


def test_log_likelihood_sanity(posterior):
    """Test log likelihood values are finite."""
    results = {}

    if 'log_likelihood' not in posterior.dtype.names:
        return results

    log_like = posterior['log_likelihood']
    n_finite = int(np.sum(np.isfinite(log_like)))
    n_total = len(log_like)

    status = "PASS" if n_finite == n_total else "FAIL"

    finite_vals = log_like[np.isfinite(log_like)]
    results['log_likelihood'] = {
        'status': status,
        'n_non_finite': n_total - n_finite,
        'max': float(np.max(finite_vals)) if n_finite > 0 else None,
        'mean': float(np.mean(finite_vals)) if n_finite > 0 else None,
        'std': float(np.std(finite_vals)) if n_finite > 0 else None,
    }

    return results


def test_railing(posterior, nbin=RAILING_NBIN, tolerance=RAILING_TOLERANCE):
    """Check for railing in luminosity_distance (upper) and chirp_mass (upper+lower)."""
    checks = {
        'luminosity_distance': ('upper',),
        'comoving_distance': ('upper', 'lower'),
        'chirp_mass': ('upper', 'lower'),
    }

    warnings = []

    for param, directions in checks.items():
        if param not in posterior.dtype.names:
            continue

        railing_down, railing_up, _, bounds = run_analyzer.railing_check(
            posterior[param], Nbin=nbin, tolerance=tolerance, label=param
        )

        if 'upper' in directions and railing_up:
            warnings.append(f"{param}: railing at upper bound (recommended max: {bounds[1]})")
        if 'lower' in directions and railing_down:
            warnings.append(f"{param}: railing at lower bound (recommended min: {bounds[0]})")

    return {
        'status': "WARNING" if warnings else "PASS",
        'warnings': warnings,
    }


def test_snr_negative(posterior):
    """Check that no SNR parameter contains negative samples."""
    snr_params = [
        name for name in posterior.dtype.names
        if ('matched_filter_snr' in name.lower() or 'optimal_snr' in name.lower()) and 'angle' not in name.lower()
    ]

    if not snr_params:
        return {}

    offenders = {}
    for param in snr_params:
        data = posterior[param]
        n_neg = int(np.sum(data < 0))
        if n_neg > 0:
            offenders[param] = {
                'n_negative': n_neg,
                'min_value': float(np.min(data)),
            }

    return {
        'status': 'WARNING' if offenders else 'PASS',
        'offenders': offenders,
    }


def test_network_snr(posterior):
    """Test network SNR consistency."""
    if 'network_optimal_snr' not in posterior.dtype.names:
        return {}

    network_snr = posterior['network_optimal_snr']

    detector_snrs = {}
    for det in ['H1', 'L1', 'V1', 'K1']:
        key = f'{det}_optimal_snr'
        if key in posterior.dtype.names:
            detector_snrs[det] = posterior[key]

    if not detector_snrs:
        return {}

    expected = np.sqrt(sum(snr ** 2 for snr in detector_snrs.values()))
    rel_error = np.abs(network_snr - expected) / (network_snr + 1e-10)
    max_error = float(np.max(rel_error))
    status = "PASS" if max_error < 0.01 else "FAIL"

    return {
        'status': status,
        'max_error': max_error,
        'network_snr_mean': float(np.mean(network_snr)),
        'detector_snrs': {det: float(np.mean(snr)) for det, snr in detector_snrs.items()},
    }


# ---------------------------------------------------------------------------
# Collapse raw test output → {status, notes}
# ---------------------------------------------------------------------------

_STATUS_RANK = {"ERROR": 4, "FAIL": 3, "WARNING": 2, "INTERESTING": 1, "PASS": 0}


def _worst(*statuses):
    """Return the most severe status from a list."""
    return max(statuses, key=lambda s: _STATUS_RANK.get(s, -1), default="N/A")


def _collapse_spin_prior(raw):
    """Collapse spin-prior test results."""
    if not raw:
        return {"status": "N/A", "notes": "no prior available"}
    statuses = [v['status'] for v in raw.values()]
    notes = []
    for v in raw.values():
        notes.extend(v.get('warnings', []))
    return {"status": _worst(*statuses), "notes": "; ".join(notes)}


def _collapse_sample_count(raw):
    """Collapse sample-count result."""
    status = raw.get('status', 'N/A')
    n = raw.get('n_samples', '?')
    notes = ""
    if status == "FAIL":
        notes = f"only {n} samples (< {raw.get('fail_threshold', 2000)} required)"
    elif status == "WARNING":
        notes = f"{n} samples (< {raw.get('warning_threshold', 5000)} recommended)"
    return {"status": status, "notes": notes}


def _collapse_multimodality(raw):
    """Collapse multimodality results."""
    if not raw:
        return {"status": "N/A", "notes": "no relevant parameters found"}
    statuses = [v['status'] for v in raw.values()]
    notes = []
    for param, v in raw.items():
        if v['status'] == 'WARNING':
            frac = v.get('max_gap_fraction', 0)
            notes.append(f"{param} gap={frac*100:.1f}% of range")
    return {"status": _worst(*statuses), "notes": "; ".join(notes)}


def _collapse_param_consistency(raw):
    """Collapse parameter-consistency results."""
    if not raw:
        return {"status": "N/A", "notes": "no derived parameters found"}
    statuses = [v['status'] for v in raw.values()]
    notes = []
    for param, v in raw.items():
        if v['status'] == 'FAIL':
            notes.append(f"{param} max_err={v.get('max_error', '?'):.2e}")
    return {"status": _worst(*statuses), "notes": "; ".join(notes)}


def _collapse_log_likelihood(raw):
    """Collapse log-evidence results."""
    if not raw:
        return {"status": "N/A", "notes": "log_likelihood not present"}
    v = raw.get('log_likelihood', {})
    status = v.get('status', 'N/A')
    n_bad = v.get('n_non_finite', 0)
    notes = f"{n_bad} non-finite log_likelihood values" if n_bad else ""
    return {"status": status, "notes": notes}


def _collapse_network_snr(raw):
    """Collapse network-SNR results."""
    if not raw:
        return {"status": "N/A", "notes": "network_optimal_snr not present"}
    status = raw.get('status', 'N/A')
    err = raw.get('max_error', None)
    notes = f"max rel. error {err*100:.2f}%" if (status == 'FAIL' and err is not None) else ""
    return {"status": status, "notes": notes}


def _collapse_diptest_unimodality(raw):
    """Collapse dip-test unimodality results."""
    if not raw:
        return {"status": "N/A", "notes": "no relevant parameters found"}
    statuses = [v['status'] for v in raw.values()]
    notes = []
    for param, v in raw.items():
        if v['status'] == 'WARNING' and v['pval'] is not None:
            notes.append(f"{param} dip={v['dip']:.4f} p={v['pval']:.3f}")
    return {"status": _worst(*statuses), "notes": "; ".join(notes)}


def _collapse_railing(raw):
    """Collapse railing results."""
    if not raw:
        return {"status": "N/A", "notes": "no relevant parameters found"}
    status = raw.get('status', 'N/A')
    notes = "; ".join(raw.get('warnings', []))
    return {"status": status, "notes": notes}


def _collapse_snr_negative(raw):
    """Collapse SNR negativity results."""
    if not raw:
        return {"status": "N/A", "notes": "no SNR parameters found"}
    status = raw.get('status', 'N/A')
    offenders = raw.get('offenders', {})
    notes = "; ".join(
        f"{p}: {v['n_negative']} negative (min={v['min_value']:.3g})"
        for p, v in offenders.items()
    )
    return {"status": status, "notes": notes}


_COLLAPSE = {
    'spin_prior': _collapse_spin_prior,
    'sample_count': _collapse_sample_count,
    'multimodality': _collapse_multimodality,
    'diptest_unimodality': _collapse_diptest_unimodality,
    'param_consistency': _collapse_param_consistency,
    'log_likelihood': _collapse_log_likelihood,
    'network_snr': _collapse_network_snr,
    'railing': _collapse_railing,
    'snr_negative': _collapse_snr_negative,
}

_TEST_LABELS = {
    'spin_prior': 'Spin Prior',
    'sample_count': 'Samples',
    'multimodality': 'Multimodality',
    'diptest_unimodality': 'Dip Test',
    'param_consistency': 'Param Consistency',
    'log_likelihood': 'Log Likelihood',
    'network_snr': 'Network SNR',
    'railing': 'Railing',
    'snr_negative': 'SNR Negative',
}


# ---------------------------------------------------------------------------
# Run all tests on one file
# ---------------------------------------------------------------------------

def run_tests_on_file(h5file, verbose=True,
                      railing_nbin=RAILING_NBIN, railing_tolerance=RAILING_TOLERANCE):
    """
    Run all quality checks on *h5file* and return a structured result dict.

    Returns
    -------
    dict with keys:
        file, timestamp, n_samples, load_error,
        tests: {test_name: {status, notes}, ...},
        overall_status
    """
    record = {
        "file": str(h5file),
        "timestamp": datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        "n_samples": None,
        "load_error": None,
        "tests": {},
        "overall_status": "N/A",
    }

    if verbose:
        print(f"\n{'='*60}")
        print(f"  {h5file}")
        print(f"{'='*60}")

    try:
        posterior, prior = open_pesummary_samples(h5file)
        record["n_samples"] = len(posterior)
        if verbose:
            n_params = len(posterior.dtype.names) if posterior.dtype.names else 0
            print(f"  Loaded {len(posterior)} samples, {n_params} parameters")
    except Exception as exc:
        record["load_error"] = str(exc)
        if verbose:
            print(f"  ERROR loading file: {exc}")
        return record

    raw_tests = {
        'spin_prior': lambda: test_spin_prior_comparison(posterior, prior),
        'sample_count': lambda: test_sample_count(posterior),
        'multimodality': lambda: test_multimodality(posterior),
        'diptest_unimodality': lambda: test_diptest_unimodality(posterior),
        'param_consistency': lambda: test_parameter_consistency(posterior),
        'log_likelihood': lambda: test_log_likelihood_sanity(posterior),
        'network_snr': lambda: test_network_snr(posterior),
        'railing': lambda: test_railing(posterior, nbin=railing_nbin, tolerance=railing_tolerance),
        'snr_negative': lambda: test_snr_negative(posterior),
    }

    for name, fn in raw_tests.items():
        try:
            raw = fn()
            collapsed = _COLLAPSE[name](raw)
        except Exception as exc:
            collapsed = {"status": "ERROR", "notes": str(exc)}
        record["tests"][name] = collapsed
        if verbose:
            flag = collapsed["status"]
            detail = f" — {collapsed['notes']}" if collapsed["notes"] else ""
            print(f"  [{flag:^7}] {_TEST_LABELS[name]}{detail}")

    statuses = [v["status"] for v in record["tests"].values()]
    record["overall_status"] = _worst(*statuses)

    if verbose:
        print(f"  Overall: {record['overall_status']}")

    return record


# ---------------------------------------------------------------------------
# Discovery helpers
# ---------------------------------------------------------------------------

_DEFAULT_EXCLUDE = {"generate-psd"}


def discover_files(root, analysis_patterns=None, exclude=None):
    """
    Walk *root* looking for posterior_samples.h5 files.

    Expected layout:
        root / {event_id} / {analysis} / posterior_samples.h5

    Parameters
    ----------
    root : str or Path
        Root results directory; each sub-directory is treated as an event.
    analysis_patterns : list[str] or None
        Names or glob patterns to match against analysis sub-directory names.
        When ``None`` or empty all analysis directories are included (subject
        to *exclude*).  Each entry may be an exact name or a glob pattern
        (e.g. ``"bilby-*"``).
    exclude : set[str] or None
        Analysis directory names to always skip.  Defaults to
        ``{"generate-psd"}``.

    Yields
    ------
    (event_id, analysis, path)
    """
    root = Path(root)
    if not root.is_dir():
        raise FileNotFoundError(f"Root directory not found: {root}")

    if exclude is None:
        exclude = _DEFAULT_EXCLUDE

    for event_dir in sorted(root.iterdir()):
        if not event_dir.is_dir():
            continue
        event_id = event_dir.name

        for analysis_dir in sorted(event_dir.iterdir()):
            if not analysis_dir.is_dir():
                continue
            analysis = analysis_dir.name

            if analysis in exclude:
                continue

            # If patterns were supplied, the analysis must match at least one.
            if analysis_patterns:
                if not any(fnmatch.fnmatch(analysis, p) for p in analysis_patterns):
                    continue

            h5file = analysis_dir / "posterior_samples.h5"
            if h5file.exists():
                yield event_id, analysis, h5file


def parse_path_ids(h5file):
    """
    Extract (event_id, analysis) from a path like:
        .../results/{event_id}/{analysis}/posterior_samples.h5
    Falls back to the parent directory name if the layout doesn't match.
    """
    p = Path(h5file).resolve()
    analysis = p.parent.name
    event_id = p.parent.parent.name
    return event_id, analysis


# ---------------------------------------------------------------------------
# Markdown report generation
# ---------------------------------------------------------------------------

_STATUS_TEXT = {
    "PASS": "PASS",
    "INTERESTING": "INTRST",
    "WARNING": "WARN",
    "FAIL": "FAIL",
    "ERROR": "ERR",
    "N/A": "N/A",
}

_STATUS_EMOJI = {
    "PASS": "✅",
    "INTERESTING": "🔍",
    "WARNING": "⚠️",
    "FAIL": "❌",
    "ERROR": "⛔",
    "N/A": "N/A",
}

_TEST_ORDER = [
    'sample_count',
    'railing',
    'snr_negative',
    'spin_prior',
    'multimodality',
    'diptest_unimodality',
    'param_consistency',
    'log_likelihood',
    'network_snr',
]


def generate_markdown(data, emoji=True, ignore=None):
    """
    Generate a markdown summary table from the JSON data produced by 'run'.

    Parameters
    ----------
    data : dict
        JSON data as returned by the run command.
    emoji : bool
        Use emoji symbols for status (default True). Pass False for plain text.
    ignore : list[str] or None
        Test names to omit from the report columns. All tests are still present
        in the underlying JSON; this only affects what is rendered.

    Returns the markdown string.
    """
    status_map = _STATUS_EMOJI if emoji else _STATUS_TEXT

    lines = []
    meta = data.get("meta", {})
    events = data.get("events", {})

    lines.append("# Automated PE Review")
    if meta.get("generated"):
        lines.append(f"\nGenerated: {meta['generated']}")
    if meta.get("root"):
        lines.append(f"Root: `{meta['root']}`")
    patterns = meta.get("analysis_patterns") or meta.get("analysis_pattern")
    if patterns:
        pat_str = ", ".join(f"`{p}`" for p in (patterns if isinstance(patterns, list) else [patterns]))
        lines.append(f"Analysis patterns: {pat_str}")
    excl = meta.get("exclude")
    if excl:
        lines.append(f"Excluded dirs: {', '.join(f'`{e}`' for e in excl)}")
    lines.append("")

    # Build flat list of rows: (event_id, analysis, record)
    rows = []
    for event_id, analyses in sorted(events.items()):
        for analysis, record in sorted(analyses.items()):
            rows.append((event_id, analysis, record))

    if not rows:
        lines.append("_No results found._")
        return "\n".join(lines)

    # Railing summary — list every event+analysis that has a non-PASS railing result
    railing_rows = [
        (event_id, analysis, record)
        for event_id, analysis, record in rows
        if not record.get("load_error")
        and record.get("tests", {}).get("railing", {}).get("status", "N/A")
        not in ("PASS", "N/A")
    ]
    lines.append("## Railing summary")
    if railing_rows:
        lines.append("")
        lines.append("The following analyses have railing warnings — review prior bounds:")
        lines.append("")
        for event_id, analysis, record in railing_rows:
            railing_notes = record["tests"]["railing"].get("notes", "")
            note_str = f" — {railing_notes}" if railing_notes else ""
            lines.append(f"- **{event_id}** / `{analysis}`{note_str}")
    else:
        lines.append("")
        lines.append("_No railing issues detected._")
    lines.append("")

    # Table header — honour --ignore (render only, JSON is unaffected)
    ignore_set = set(ignore or [])
    test_cols = [t for t in _TEST_ORDER if t not in ignore_set]
    header_parts = ["Event", "Analysis", "N samples", "Overall"]
    header_parts += [_TEST_LABELS[t] for t in test_cols]
    header_parts.append("Notes")

    sep_parts = ["---"] * len(header_parts)

    lines.append("| " + " | ".join(header_parts) + " |")
    lines.append("| " + " | ".join(sep_parts) + " |")

    for event_id, analysis, record in rows:
        if record.get("load_error"):
            row = [event_id, analysis, "—", status_map.get("ERROR", "ERR")]
            row += ["—"] * len(test_cols)
            row.append(f"Load error: {record['load_error']}")
            lines.append("| " + " | ".join(row) + " |")
            continue

        n_samples = record.get("n_samples")
        n_str = f"{n_samples:,}" if n_samples is not None else "—"
        overall = record.get("overall_status", "N/A")

        row = [event_id, analysis, n_str, status_map.get(overall, overall)]

        # Collect notes across all tests
        all_notes = []
        for t in test_cols:
            test_result = record.get("tests", {}).get(t, {})
            status = test_result.get("status", "N/A")
            row.append(status_map.get(status, status))
            note = test_result.get("notes", "")
            if note:
                all_notes.append(note)

        row.append("; ".join(all_notes) if all_notes else "")
        lines.append("| " + " | ".join(row) + " |")

    # Summary counts
    n_pass = sum(1 for _, _, r in rows if r.get("overall_status") == "PASS")
    n_interesting = sum(1 for _, _, r in rows if r.get("overall_status") == "INTERESTING")
    n_warn = sum(1 for _, _, r in rows if r.get("overall_status") == "WARNING")
    n_fail = sum(1 for _, _, r in rows if r.get("overall_status") == "FAIL")
    n_err = sum(1 for _, _, r in rows if r.get("load_error"))

    lines.append("")
    lines.append(f"**Total runs:** {len(rows)} | "
                 f"PASS: {n_pass} | INTERESTING: {n_interesting} | "
                 f"WARNING: {n_warn} | FAIL: {n_fail} | ERROR: {n_err}")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _collect_h5files(args):
    """Resolve the list of (event_id, analysis, path) tuples from parsed args."""
    if getattr(args, 'h5files', None):
        return [(eid, an, Path(f))
                for f in args.h5files
                for eid, an in [parse_path_ids(f)]]
    if args.root:
        patterns = args.analysis or []  # empty list → discover all
        exclude = set(args.exclude) if getattr(args, 'exclude', None) else None
        if patterns:
            print(f"Discovering analyses matching {patterns} under {args.root} ...")
        else:
            excl_note = f" (excluding: {sorted(exclude or _DEFAULT_EXCLUDE)})"
            print(f"Discovering all analyses under {args.root}{excl_note} ...")
        found = list(discover_files(args.root, analysis_patterns=patterns, exclude=exclude))
        if not found:
            print("No matching posterior_samples.h5 files found.")
            return None
        print(f"Found {len(found)} run(s).")
        return found
    print("ERROR: Supply either H5 file paths or --root.")
    return None


def _do_run(args):
    """Run tests and return the result data dict (does not write markdown)."""
    h5files = _collect_h5files(args)
    if h5files is None:
        return None, 1

    events = {}
    for event_id, analysis, path in h5files:
        record = run_tests_on_file(
            path,
            verbose=not args.quiet,
            railing_nbin=getattr(args, 'railing_nbin', RAILING_NBIN),
            railing_tolerance=getattr(args, 'railing_tolerance', RAILING_TOLERANCE),
        )
        events.setdefault(event_id, {})[analysis] = record

    data = {
        "meta": {
            "generated": datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            "root": str(args.root) if args.root else None,
            "analysis_patterns": args.analysis or [],
            "exclude": list(args.exclude) if getattr(args, 'exclude', None) else sorted(_DEFAULT_EXCLUDE),
        },
        "events": events,
    }

    out_path = args.output or "review.json"
    with open(out_path, 'w') as fh:
        json.dump(data, fh, indent=2)
    print(f"\nResults saved to {out_path}")
    return data, 0


def _do_report(data, markdown_out, emoji=True, ignore=None):
    """Write markdown from *data*. Prints to stdout if *markdown_out* is None."""
    md = generate_markdown(data, emoji=emoji, ignore=ignore)
    if markdown_out:
        with open(markdown_out, 'w') as fh:
            fh.write(md)
            fh.write("\n")
        print(f"Markdown report written to {markdown_out}")
    else:
        print(md)


def cmd_run(args):
    """'run' subcommand: run tests, write JSON, optionally write markdown."""
    data, rc = _do_run(args)
    if rc != 0:
        return rc
    if getattr(args, 'markdown', None):
        _do_report(data, args.markdown, emoji=not args.no_emoji,
                   ignore=getattr(args, 'ignore', None))
    return 0


def cmd_report(args):
    """'report' subcommand: generate markdown from an existing JSON file."""
    with open(args.json_file) as fh:
        data = json.load(fh)
    _do_report(data, args.output, emoji=not args.no_emoji,
               ignore=getattr(args, 'ignore', None))
    return 0


def cmd_run_and_report(args):
    """Default (no subcommand): run tests then immediately write markdown."""
    data, rc = _do_run(args)
    if rc != 0:
        return rc
    _do_report(data, args.markdown, emoji=not args.no_emoji,
               ignore=getattr(args, 'ignore', None))
    return 0


def _add_run_args(p, with_h5files=False):
    """Attach the run-mode arguments to parser *p*.

    Parameters
    ----------
    with_h5files : bool
        Include the positional H5FILE argument. Only set for the 'run'
        subcommand; the default top-level mode uses --root/--analysis only,
        which avoids ambiguity with the 'report' subcommand name.
    """
    if with_h5files:
        p.add_argument(
            'h5files',
            nargs='*',
            metavar='H5FILE',
            help='Explicit posterior_samples.h5 paths (alternative to --root).',
        )
    p.add_argument(
        '--root', '-r',
        metavar='DIR',
        help='Root results directory. Events are sub-directories of this path.',
    )
    p.add_argument(
        '--analysis', '-a',
        nargs='*',
        metavar='PATTERN',
        default=None,
        help='Analysis names or glob patterns to include (default: all). '
             'Pass one or more values, e.g. '
             '"bilby-IMRPhenomXPHM-SpinTaylor bilby-IMRPhenomXPNR" or "bilby-*". '
             'Omit entirely to include every analysis directory.',
    )
    p.add_argument(
        '--exclude',
        nargs='*',
        metavar='DIR',
        default=sorted(_DEFAULT_EXCLUDE),
        help='Analysis directory names to always skip '
             f'(default: {sorted(_DEFAULT_EXCLUDE)}). '
             'Pass an empty list (--exclude) to disable all exclusions.',
    )
    p.add_argument(
        '--output', '-o',
        metavar='JSON',
        default='review.json',
        help='Output JSON file (default: review.json).',
    )
    p.add_argument(
        '--markdown', '-m',
        metavar='MD',
        help='Also write a markdown report to this file '
             '(default: print to stdout when no subcommand, skip when using run).',
    )
    p.add_argument(
        '--no-emoji',
        action='store_true',
        default=False,
        help='Use plain text (PASS/WARN/FAIL) instead of emoji in the markdown report.',
    )
    p.add_argument(
        '--ignore',
        nargs='+',
        metavar='TEST',
        default=None,
        help='Test names to omit from the markdown report (all tests still run '
             'and are saved to JSON). E.g. --ignore multimodality spin_prior. '
             f'Valid names: {", ".join(_TEST_ORDER)}.',
    )
    p.add_argument(
        '--railing-nbin',
        type=int,
        default=RAILING_NBIN,
        metavar='N',
        help=f'Number of bins for the railing check (default: {RAILING_NBIN}).',
    )
    p.add_argument(
        '--railing-tolerance',
        type=float,
        default=RAILING_TOLERANCE,
        metavar='TOL',
        help=f'Tolerance for the railing check (default: {RAILING_TOLERANCE}).',
    )
    p.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Suppress per-test console output.',
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            'Batch PE review. '
            'With no subcommand: runs tests and prints a markdown report. '
            'Use "run" to only produce JSON, "report" to render an existing JSON.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Default (no subcommand): --root/--analysis only, no positional h5files.
    # Positional args at the top level conflict with subcommand name detection.
    _add_run_args(parser, with_h5files=False)

    sub = parser.add_subparsers(dest='command')

    # --- run subcommand ---
    p_run = sub.add_parser(
        'run',
        help='Run quality checks and write JSON (use --markdown to also write a report).',
    )
    _add_run_args(p_run, with_h5files=True)

    # --- report subcommand ---
    p_rep = sub.add_parser(
        'report',
        help='Generate a markdown table from an existing JSON results file.',
    )
    p_rep.add_argument(
        'json_file',
        metavar='JSON',
        help='JSON file produced by the run subcommand.',
    )
    p_rep.add_argument(
        '--output', '-o',
        metavar='MD',
        help='Output markdown file (default: print to stdout).',
    )
    p_rep.add_argument(
        '--no-emoji',
        action='store_true',
        default=False,
        help='Use plain text (PASS/WARN/FAIL) instead of emoji.',
    )
    p_rep.add_argument(
        '--ignore',
        nargs='+',
        metavar='TEST',
        default=None,
        help='Test names to omit from the markdown report (data in JSON is unaffected). '
             f'Valid names: {", ".join(_TEST_ORDER)}.',
    )

    args = parser.parse_args()

    if args.command == 'run':
        return cmd_run(args)
    if args.command == 'report':
        return cmd_report(args)
    # No subcommand: run + report
    return cmd_run_and_report(args)


if __name__ == '__main__':
    sys.exit(main())
