#!/usr/bin/env python3
"""Verify that all expected release files exist and are readable.

For each event in the ledger that has at least one APPROVED bilby production,
checks the Stage 2B combined release:
  1. releases/{release}/{event}/C00-Combined.h5 exists and is non-empty
  2. Every expected label (C00:{waveform}) is present in the HDF5
  3. Each present label has at least one posterior sample
  4. Each expected skymap FITS file exists and is non-empty

Also optionally checks the Stage 1 SpinTaylor-only release:
  releases/{release}/spintaylor/{event}/C00-Combined.h5

Usage:
    python tests/test_release_files.py [--ledger .asimov/ledger.yml]
                                       [--release-dir releases/staging]
                                       [--spintaylor-dir releases/staging/spintaylor]
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from ledger import (
    load_ledger, get_events, approved_bilby_productions, release_label,
)


import h5py


def _expected_labels(productions):
    labels = [
        release_label(settings["waveform"]["approximant"])
        for _, settings in productions
    ]
    if len(labels) > 1:
        labels.append("Mixed")
    return labels


def check_event(event_name, event, release_dir):
    """Yield (description, passed, detail) tuples for one event."""
    productions = list(approved_bilby_productions(event))
    if not productions:
        return

    expected_labels = _expected_labels(productions)
    hdf5_path = os.path.join(release_dir, event_name, "C00-Combined.h5")

    exists = os.path.isfile(hdf5_path)
    yield (f"{hdf5_path} exists", exists, "")
    if not exists:
        return

    size = os.path.getsize(hdf5_path)
    yield (f"{hdf5_path} non-empty", size > 0, f"{size} bytes")
    if not size:
        return

    try:
        with h5py.File(hdf5_path, "r") as fh:
            for label in expected_labels:
                key = f"C00:{label}"
                found = key in fh
                present = sorted(k for k in fh.keys() if k not in ("history", "version"))
                yield (
                    f"{event_name}: label '{key}' present",
                    found,
                    "" if found else f"found: {present}",
                )
                if found:
                    try:
                        n = fh[key]["posterior_samples"].shape[0]
                        yield (f"{event_name}: '{key}' has samples", n > 0, f"{n} samples")
                    except KeyError as exc:
                        yield (f"{event_name}: '{key}' posterior_samples readable", False, str(exc))
    except Exception as exc:
        yield (f"{hdf5_path} readable", False, str(exc))
        return

    for label in expected_labels:
        fits_name = f"C00-{label.replace(':', '-')}.fits.gz"
        fits_path = os.path.join(release_dir, event_name, fits_name)
        fits_ok = os.path.isfile(fits_path) and os.path.getsize(fits_path) > 0
        detail = f"{os.path.getsize(fits_path)} bytes" if os.path.isfile(fits_path) else "missing"
        yield (f"{event_name}: skymap {fits_name} exists", fits_ok, detail)


def check_spintaylor_event(event_name, event, spintaylor_dir):
    """Yield (description, passed, detail) tuples for one Stage 1 SpinTaylor event."""
    # Only events with exactly one approved IMRPhenomXPHM production are in Stage 1.
    qualifying = [
        production
        for production, settings in approved_bilby_productions(event)
        if settings["waveform"]["approximant"] == "IMRPhenomXPHM"
    ]
    if len(qualifying) != 1:
        return

    hdf5_path = os.path.join(spintaylor_dir, event_name, "C00-Combined.h5")
    exists = os.path.isfile(hdf5_path)
    yield (f"{hdf5_path} exists", exists, "")
    if not exists:
        return

    size = os.path.getsize(hdf5_path)
    yield (f"{hdf5_path} non-empty", size > 0, f"{size} bytes")
    if not size:
        return

    try:
        with h5py.File(hdf5_path, "r") as fh:
            key = "C00:IMRPhenomXPHM-SpinTaylor"
            found = key in fh
            present = sorted(k for k in fh.keys() if k not in ("history", "version"))
            yield (
                f"{event_name}: label '{key}' present",
                found,
                "" if found else f"found: {present}",
            )
            if found:
                try:
                    n = fh[key]["posterior_samples"].shape[0]
                    yield (f"{event_name}: '{key}' has samples", n > 0, f"{n} samples")
                except KeyError as exc:
                    yield (f"{event_name}: '{key}' posterior_samples readable", False, str(exc))
    except Exception as exc:
        yield (f"{hdf5_path} readable", False, str(exc))


def main():
    parser = argparse.ArgumentParser(description="Verify release file completeness and readability")
    parser.add_argument("--ledger", default=".asimov/ledger.yml")
    parser.add_argument("--release-dir", default="releases/staging")
    parser.add_argument(
        "--spintaylor-dir",
        default=None,
        help="Stage 1 SpinTaylor release directory (default: {release-dir}/spintaylor if it exists)",
    )
    args = parser.parse_args()

    spintaylor_dir = args.spintaylor_dir
    if spintaylor_dir is None:
        candidate = os.path.join(args.release_dir, "spintaylor")
        if os.path.isdir(candidate):
            spintaylor_dir = candidate

    data = load_ledger(args.ledger)
    events = get_events(data)

    passed = failed = 0
    print(f"Stage 2B combined release checks ({args.release_dir}):")
    for name, event in sorted(events.items()):
        for description, ok, detail in check_event(name, event, args.release_dir):
            status = "PASS" if ok else "FAIL"
            suffix = f"  ({detail})" if detail else ""
            print(f"  [{status}] {description}{suffix}")
            if ok:
                passed += 1
            else:
                failed += 1

    if spintaylor_dir:
        print(f"\nStage 1 SpinTaylor checks ({spintaylor_dir}):")
        for name, event in sorted(events.items()):
            for description, ok, detail in check_spintaylor_event(name, event, spintaylor_dir):
                status = "PASS" if ok else "FAIL"
                suffix = f"  ({detail})" if detail else ""
                print(f"  [{status}] {description}{suffix}")
                if ok:
                    passed += 1
                else:
                    failed += 1

    total = passed + failed
    print(f"\n{total} checks: {passed} passed, {failed} failed")
    sys.exit(0 if failed == 0 else 1)


if __name__ == "__main__":
    main()


# ============================================================================
# pytest interface
# ============================================================================
# The functions above provide a standalone CLI.  The functions below expose
# the same checks as individual pytest test functions so that
#
#     pytest tests/test_release_files.py \
#         --ledger project/.asimov/ledger.yml \
#         --release-dir releases/v4.0
#
# gives per-event, per-check output.  The event_name_and_data fixture is
# injected by conftest.pytest_generate_tests.

try:
    import pytest as _pytest

    def test_hdf5_present(event_name_and_data, release_dir):
        """C00-Combined.h5 must exist and be non-empty."""
        name, event = event_name_and_data
        hdf5 = os.path.join(release_dir, name, "C00-Combined.h5")
        assert os.path.isfile(hdf5), f"Missing: {hdf5}"
        assert os.path.getsize(hdf5) > 0, f"Empty: {hdf5}"

    def test_hdf5_labels_and_samples(event_name_and_data, release_dir):
        """All expected waveform labels must be present with posterior samples."""
        name, event = event_name_and_data
        productions = list(approved_bilby_productions(event))
        expected = _expected_labels(productions)
        hdf5 = os.path.join(release_dir, name, "C00-Combined.h5")
        if not os.path.isfile(hdf5):
            _pytest.skip(f"HDF5 not found: {hdf5}")
        failures = []
        with h5py.File(hdf5, "r") as fh:
            for label in expected:
                key = f"C00:{label}"
                if key not in fh:
                    failures.append(f"{key!r} missing")
                    continue
                try:
                    n = fh[key]["posterior_samples"].shape[0]
                    if n == 0:
                        failures.append(f"{key!r}: 0 samples")
                except KeyError:
                    failures.append(f"{key!r}: posterior_samples unreadable")
        assert not failures, "Label/sample failures: " + "; ".join(failures)

    def test_skymaps_present(event_name_and_data, release_dir):
        """Per-label skymap FITS files must exist and be non-empty."""
        name, event = event_name_and_data
        productions = list(approved_bilby_productions(event))
        expected = _expected_labels(productions)
        missing = []
        for label in expected:
            fits_name = f"C00-{label.replace(':', '-')}.fits.gz"
            fits_path = os.path.join(release_dir, name, fits_name)
            if not os.path.isfile(fits_path) or os.path.getsize(fits_path) == 0:
                missing.append(fits_name)
        assert not missing, "Missing/empty skymaps: " + ", ".join(missing)

except ImportError:
    pass
