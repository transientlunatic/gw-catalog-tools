"""Quality checks for release PE files.

Tests run against the populated release directory:

    pytest tests/test_release_quality.py \\
        --ledger project/.asimov/ledger.yml \\
        --release-dir releases/v4.0

Checks:
  - Physical parameter bounds (masses, spins, distance, mass ratio)
  - Geocentric time consistent with event trigger time
  - Skymap FITS files readable with normalised probability
  - Cross-waveform JS divergence on key intrinsic parameters
  - C00:Mixed samples lie within the union of the per-waveform posteriors

The ``event_name_and_data`` fixture is parametrized per event by
conftest.pytest_generate_tests.
"""

import os
import sys
import warnings

import h5py
import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from ledger import approved_bilby_productions, release_label


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _combined_hdf5(release_dir, event_name):
    return os.path.join(release_dir, event_name, "C00-Combined.h5")


def _waveform_labels(event):
    """Return the C00-prefixed waveform labels for approved productions."""
    return [
        f"C00:{release_label(s['waveform']['approximant'])}"
        for _, s in approved_bilby_productions(event)
    ]


def _open_samples(fh, key):
    return np.array(fh[key]["posterior_samples"])


def _js_divergence(a, b, bins=100):
    """Jensen-Shannon divergence in bits between two 1-D sample arrays."""
    lo = min(a.min(), b.min())
    hi = max(a.max(), b.max())
    if lo == hi:
        return 0.0
    p, _ = np.histogram(a, bins=bins, range=(lo, hi))
    q, _ = np.histogram(b, bins=bins, range=(lo, hi))
    # Laplace smoothing to avoid log(0)
    p = (p + 1).astype(float)
    q = (q + 1).astype(float)
    p /= p.sum()
    q /= q.sum()
    m = 0.5 * (p + q)
    return 0.5 * (np.sum(p * np.log2(p / m)) + np.sum(q * np.log2(q / m)))


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_physical_bounds(event_name_and_data, release_dir):
    """All posterior samples must lie within physical parameter bounds.

    Checks per waveform label in the combined HDF5:
      - mass_1, mass_2 > 0
      - mass_2 <= mass_1  (i.e. q <= 1)
      - mass_ratio in (0, 1]
      - a_1, a_2 in [0, 1]
      - luminosity_distance > 0
    """
    name, event = event_name_and_data
    path = _combined_hdf5(release_dir, name)
    if not os.path.isfile(path):
        pytest.skip(f"HDF5 not found: {path}")

    violations = []
    with h5py.File(path, "r") as fh:
        analysis_keys = [k for k in fh if k not in ("history", "version")]
        for key in analysis_keys:
            s = _open_samples(fh, key)
            cols = s.dtype.names or []

            for mass in ("mass_1", "mass_2"):
                if mass in cols and np.any(s[mass] <= 0):
                    violations.append(f"{key}/{mass}: non-positive values")

            if "mass_1" in cols and "mass_2" in cols:
                if np.any(s["mass_2"] > s["mass_1"]):
                    violations.append(f"{key}: mass_2 > mass_1 (mass ratio > 1)")

            if "mass_ratio" in cols:
                if np.any(s["mass_ratio"] <= 0) or np.any(s["mass_ratio"] > 1):
                    violations.append(f"{key}/mass_ratio: values outside (0, 1]")

            for spin in ("a_1", "a_2"):
                if spin in cols and (np.any(s[spin] < 0) or np.any(s[spin] > 1)):
                    violations.append(f"{key}/{spin}: outside [0, 1]")

            if "luminosity_distance" in cols and np.any(s["luminosity_distance"] <= 0):
                violations.append(f"{key}/luminosity_distance: non-positive values")

    assert not violations, "Physical bound violations:\n" + "\n".join(violations)


def test_geocentric_time(event_name_and_data, release_dir):
    """Posterior median geocent_time must be within 0.1 s of the trigger time.

    Skips if the ledger does not record a trigger time for the event.
    """
    name, event = event_name_and_data
    trigger_time = (
        event.get("trigger_time")
        or event.get("gps_time")
        or event.get("trigger-time")
    )
    if trigger_time is None:
        pytest.skip("No trigger_time in ledger for this event")

    path = _combined_hdf5(release_dir, name)
    if not os.path.isfile(path):
        pytest.skip(f"HDF5 not found: {path}")

    violations = []
    with h5py.File(path, "r") as fh:
        for key in [k for k in fh if k not in ("history", "version")]:
            s = _open_samples(fh, key)
            if "geocent_time" not in (s.dtype.names or []):
                continue
            median_t = np.median(s["geocent_time"])
            offset = abs(median_t - float(trigger_time))
            if offset > 0.1:
                violations.append(
                    f"{key}: median geocent_time offset {offset:.3f} s "
                    f"(trigger={trigger_time}, median={median_t:.3f})"
                )

    assert not violations, "geocent_time sanity failures:\n" + "\n".join(violations)


def test_skymap_fits_valid(event_name_and_data, release_dir):
    """Skymap FITS files must be readable and contain a normalised probability column.

    Requires astropy; skips if not installed.
    """
    astropy_fits = pytest.importorskip("astropy.io.fits",
                                       reason="astropy required for FITS validation")
    name, event = event_name_and_data
    labels = [
        release_label(s["waveform"]["approximant"])
        for _, s in approved_bilby_productions(event)
    ]
    if len(labels) > 1:
        labels.append("Mixed")

    failures = []
    for label in labels:
        fits_name = f"C00-{label.replace(':', '-')}.fits.gz"
        fits_path = os.path.join(release_dir, name, fits_name)
        if not os.path.isfile(fits_path):
            failures.append(f"{fits_name}: missing")
            continue
        try:
            with astropy_fits.open(fits_path) as hdul:
                prob_cols = [
                    c.name for c in hdul[1].columns
                    if c.name.upper() in ("PROB", "PROBABILITY")
                ]
                if not prob_cols:
                    failures.append(f"{fits_name}: no PROB column")
                    continue
                total = float(hdul[1].data[prob_cols[0]].sum())
                if not 0.99 < total < 1.01:
                    failures.append(
                        f"{fits_name}: PROB sums to {total:.5f}, expected ~1.0"
                    )
        except Exception as exc:
            failures.append(f"{fits_name}: unreadable ({exc})")

    assert not failures, "Skymap FITS failures:\n" + "\n".join(failures)


def test_cross_waveform_js_divergence(event_name_and_data, release_dir):
    """JS divergence between waveform posteriors must be < 0.6 bits for key params.

    Events with only one approved waveform are skipped.

    Parameters checked: chirp_mass, mass_ratio, chi_eff, luminosity_distance.

    Threshold:
      > 0.3 bits  — emits a UserWarning (visible in pytest -W output)
      > 0.6 bits  — hard failure
    """
    name, event = event_name_and_data
    labels = _waveform_labels(event)
    if len(labels) < 2:
        pytest.skip("Only one approved waveform — no cross-waveform comparison")

    path = _combined_hdf5(release_dir, name)
    if not os.path.isfile(path):
        pytest.skip(f"HDF5 not found: {path}")

    WARN_THRESHOLD = 0.3   # bits
    FAIL_THRESHOLD = 0.6   # bits
    params = ["chirp_mass", "mass_ratio", "chi_eff", "luminosity_distance"]

    with h5py.File(path, "r") as fh:
        present = [l for l in labels if l in fh]
        if len(present) < 2:
            pytest.skip("Fewer than 2 waveform labels found in HDF5")
        samples = {l: _open_samples(fh, l) for l in present}

    warn_msgs = []
    fail_msgs = []
    for i, la in enumerate(present):
        for lb in present[i + 1:]:
            sa, sb = samples[la], samples[lb]
            for p in params:
                if p not in (sa.dtype.names or []) or p not in (sb.dtype.names or []):
                    continue
                jsd = _js_divergence(sa[p], sb[p])
                msg = f"{la} vs {lb} | {p}: JS = {jsd:.3f} bits"
                if jsd > FAIL_THRESHOLD:
                    fail_msgs.append(msg)
                elif jsd > WARN_THRESHOLD:
                    warn_msgs.append(msg)

    if warn_msgs:
        warnings.warn(
            f"{name}: cross-waveform JS divergence above 0.3 bits:\n"
            + "\n".join(warn_msgs),
            UserWarning,
            stacklevel=2,
        )

    assert not fail_msgs, (
        f"Cross-waveform JS divergence > {FAIL_THRESHOLD} bits:\n"
        + "\n".join(fail_msgs)
    )


def test_mixed_posterior_consistency(event_name_and_data, release_dir):
    """C00:Mixed samples must lie within the union of the individual waveform ranges.

    At most 0.5% of samples are allowed outside the per-waveform union range
    (to allow for floating-point edge effects from the combining step).
    Events with only one approved waveform are skipped.
    """
    name, event = event_name_and_data
    labels = _waveform_labels(event)
    if len(labels) < 2:
        pytest.skip("Single waveform — no Mixed posterior expected")

    path = _combined_hdf5(release_dir, name)
    if not os.path.isfile(path):
        pytest.skip(f"HDF5 not found: {path}")

    params = ["chirp_mass", "mass_ratio", "luminosity_distance"]
    TAIL_TOLERANCE = 0.005   # 0.5%

    with h5py.File(path, "r") as fh:
        if "C00:Mixed" not in fh:
            pytest.skip("C00:Mixed not present in HDF5")
        mixed = _open_samples(fh, "C00:Mixed")
        per_waveform = {l: _open_samples(fh, l) for l in labels if l in fh}

    issues = []
    for p in params:
        if p not in (mixed.dtype.names or []):
            continue
        wf_with_param = [s for s in per_waveform.values() if p in (s.dtype.names or [])]
        if not wf_with_param:
            continue
        union_min = min(s[p].min() for s in wf_with_param)
        union_max = max(s[p].max() for s in wf_with_param)
        frac_outside = float(np.mean(
            (mixed[p] < union_min) | (mixed[p] > union_max)
        ))
        if frac_outside > TAIL_TOLERANCE:
            issues.append(
                f"{p}: {frac_outside * 100:.1f}% of Mixed samples outside "
                f"individual-waveform range [{union_min:.4g}, {union_max:.4g}]"
            )

    assert not issues, "C00:Mixed posterior out of bounds:\n" + "\n".join(issues)
