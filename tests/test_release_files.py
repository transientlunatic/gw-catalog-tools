#!/usr/bin/env python3
"""Verify that all expected release files exist and are readable.

For each event in the ledger that has at least one APPROVED bilby production,
checks:
  1. releases/{release}/{event}/C00-Combined.h5 exists and is non-empty
  2. Every expected label (C00:{waveform}) is present in the HDF5
  3. Each present label has at least one posterior sample
  4. Each expected skymap FITS file exists and is non-empty

Usage:
    python tests/test_release_files.py [--ledger .asimov/ledger.yml]
                                       [--release-dir releases/staging]
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


def main():
    parser = argparse.ArgumentParser(description="Verify release file completeness and readability")
    parser.add_argument("--ledger", default=".asimov/ledger.yml")
    parser.add_argument("--release-dir", default="releases/staging")
    args = parser.parse_args()

    data = load_ledger(args.ledger)
    events = get_events(data)

    passed = failed = 0
    for name, event in sorted(events.items()):
        for description, ok, detail in check_event(name, event, args.release_dir):
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
