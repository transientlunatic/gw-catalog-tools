#!/usr/bin/env python3
"""
Run the automated PE review checks on one posterior samples file.

This needs neither asimov nor a project.  It runs exactly the checks of the
asimov tests (``test_automated_review.py``), because both call the functions in
``review_utils``.

Usage:
    python test_automated_review_standalone.py <path_to_posterior_samples.h5>

The exit status is 1 if any check FAILs and 0 otherwise.  Warnings do not
affect it.
"""

import argparse
import sys
from datetime import datetime

from review_utils import open_pesummary_samples, read_run_info, run_all_checks



def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Run the automated PE review checks on a single HDF5 file")
    parser.add_argument("h5file", help="Path to posterior_samples.h5 file")
    args = parser.parse_args(argv)

    print("=" * 60)
    print("AUTOMATED PE REVIEW - STANDALONE")
    print("=" * 60)
    print(f"File: {args.h5file}")
    print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    try:
        posterior, prior = open_pesummary_samples(args.h5file)
        info = read_run_info(args.h5file)
    except Exception as error:
        print(f"\nERROR: Failed to load file: {error}")
        return 1
    print(f"Loaded {len(posterior)} samples with {len(posterior.dtype.names)} parameters\n")

    results = run_all_checks(posterior, prior, info)

    counts = {}
    for name, result in results.items():
        counts[result["status"]] = counts.get(result["status"], 0) + 1
        print(f"{result['status']:8} {name.replace('_', ' '):24} {result['message']}")

    print("\n" + "=" * 60)
    print("SUMMARY: " + ", ".join(f"{count} {status}" for status, count in sorted(counts.items())))
    print("=" * 60)
    return 1 if counts.get("FAIL") else 0


if __name__ == "__main__":
    sys.exit(main())
