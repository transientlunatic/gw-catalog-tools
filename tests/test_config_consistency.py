"""Check that all analyses contributing to a combined release file used
consistent data-conditioning settings.

For each combined HDF5, the config stored under each waveform label is
compared against all other labels.  Differences in the waveform approximant,
waveform-specific arguments, and job-submission settings are expected and
excluded.  Differences in the detector network, segment duration, PSD, trigger
time, or frequency range are unexpected and warrant review.

Results are emitted as UserWarnings so they appear in pytest output without
failing the suite — an unusual config difference may be deliberate.

Run with:
    pytest tests/test_config_consistency.py \\
        --ledger project/.asimov/ledger.yml \\
        --release-dir releases/v4.0
"""

import os
import sys
import warnings

import h5py
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from ledger import approved_bilby_productions, release_label


# Keys that are expected to differ between waveform analyses.
# Everything else is compared and any discrepancy is flagged.
_WAVEFORM_SPECIFIC_KEYS = frozenset({
    # Approximant / model
    "approximant",
    "waveform-arguments-dict",
    "waveform-arguments",
    "frequency-domain-source-model",
    "time-domain-source-model",
    "waveform-generator",
    "waveform-generator-class",
    # Output / job identity
    "label",
    "outdir",
    "webdir",
    # HTCondor / scheduler
    "accounting",
    "accounting-user",
    "accounting_user",
    "request-memory",
    "request-cpus",
    "request-disk",
    "scheduler-env",
    "scheduler-analysis-time",
    "extra-kwargs",
    # NR / surrogate specific
    "nr-file",
    "NR-file",
})

# These keys are the most safety-critical if they differ — they get flagged
# by name in the output so reviewers know which to prioritise.
_HIGH_PRIORITY_KEYS = frozenset({
    "trigger-time",
    "trigger_time",
    "detectors",
    "detectors-list",
    "duration",
    "sampling-frequency",
    "minimum-frequency",
    "maximum-frequency",
    "reference-frequency",
    "channel-dict",
    "psd-dict",
})


def _read_label_config(fh, key):
    """Extract the config dict stored under a PESummary label group.

    Handles two common storage formats:
      1. Parallel datasets ``config_file/keys`` and ``config_file/vals``
      2. Individual datasets named by config key under ``config_file/``

    Returns an empty dict if no config is stored.
    """
    if key not in fh or "config_file" not in fh[key]:
        return {}

    cfg = fh[key]["config_file"]

    def _decode(x):
        return x.decode() if isinstance(x, bytes) else str(x)

    # Format 1: parallel key/val arrays
    if "keys" in cfg and "vals" in cfg:
        return {_decode(k): _decode(v) for k, v in zip(cfg["keys"], cfg["vals"])}

    # Format 2: one dataset per config key
    result = {}
    for k in cfg:
        try:
            raw = cfg[k][()]
            result[k] = raw.decode() if isinstance(raw, bytes) else str(raw)
        except Exception:
            result[k] = "<unreadable>"
    return result


def _compare_configs(configs):
    """Compare config dicts across labels.

    Parameters
    ----------
    configs : dict[str, dict]
        Mapping from label name to its config dict.

    Returns
    -------
    list of dict
        Each entry describes one differing key:
        ``{"key": str, "values": {label: value}, "high_priority": bool}``
    """
    if len(configs) < 2:
        return []

    all_keys = set()
    for cfg in configs.values():
        all_keys.update(cfg.keys())
    comparable_keys = all_keys - _WAVEFORM_SPECIFIC_KEYS

    diffs = []
    for key in sorted(comparable_keys):
        values = {label: cfg.get(key, "<absent>") for label, cfg in configs.items()}
        unique = set(values.values())
        if len(unique) > 1:
            diffs.append({
                "key": key,
                "values": values,
                "high_priority": key in _HIGH_PRIORITY_KEYS,
            })
    return diffs


def _truncate(s, n=80):
    return s if len(s) <= n else s[:n] + "…"


def test_config_consistency(event_name_and_data, release_dir):
    """Config settings shared across waveform analyses must be identical.

    Differences in waveform model, job-submission settings, and output
    directories are excluded.  All other differences emit a UserWarning.
    High-priority differences (detectors, duration, trigger time, etc.)
    cause the test to fail.
    """
    name, event = event_name_and_data
    productions = list(approved_bilby_productions(event))
    if len(productions) < 2:
        pytest.skip("Only one approved waveform — nothing to compare")

    path = os.path.join(release_dir, name, "C00-Combined.h5")
    if not os.path.isfile(path):
        pytest.skip(f"HDF5 not found: {path}")

    labels = [
        f"C00:{release_label(s['waveform']['approximant'])}"
        for _, s in productions
    ]

    with h5py.File(path, "r") as fh:
        configs = {
            label: _read_label_config(fh, label)
            for label in labels
            if label in fh
        }

    if len(configs) < 2:
        pytest.skip("Fewer than 2 waveform labels found in HDF5")

    diffs = _compare_configs(configs)
    if not diffs:
        return

    high_priority = [d for d in diffs if d["high_priority"]]
    low_priority = [d for d in diffs if not d["high_priority"]]

    def _fmt_diff(d):
        vals = "; ".join(
            f"{lbl}={_truncate(v)}" for lbl, v in sorted(d["values"].items())
        )
        priority = "[HIGH PRIORITY] " if d["high_priority"] else ""
        return f"  {priority}{d['key']}: {vals}"

    all_lines = [_fmt_diff(d) for d in diffs]

    if low_priority:
        warnings.warn(
            f"{name}: config inconsistencies (non-critical):\n"
            + "\n".join(_fmt_diff(d) for d in low_priority),
            UserWarning,
            stacklevel=2,
        )

    assert not high_priority, (
        f"{name}: high-priority config inconsistencies between waveform analyses:\n"
        + "\n".join(_fmt_diff(d) for d in high_priority)
        + "\n\nAll differences:\n"
        + "\n".join(all_lines)
    )
