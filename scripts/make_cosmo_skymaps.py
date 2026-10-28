"""
Generate missing cosmologically-reweighted skymaps for the GWTC-3 old-format release.

For each *_cosmo.h5 in the release directory, reads the waveform labels
(top-level HDF5 groups, excluding 'history' and 'version') and emits a
ligo-skymap-from-samples command for any label whose cosmo_reweight skymap
does not already exist on disk.

Usage (from any directory):
    python scripts/make_cosmo_skymaps.py \\
        --release-dir /home/pe.o3/tmp/o3b-release-2024/cosmo \\
        --ligo-skymap /path/to/ligo-skymap-from-samples \\
        [--jobs 16] [--maxpts 2000] [--dry-run]
"""

import argparse
import os
import subprocess
import sys


def get_h5_labels(h5path: str) -> list[str]:
    """Return top-level group names from an HDF5 file using h5ls."""
    result = subprocess.run(
        ["h5ls", h5path],
        capture_output=True, text=True
    )
    labels = []
    for line in result.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[1] == "Group":
            name = parts[0]
            if name not in ("history", "version"):
                labels.append(name)
    return labels


def cosmo_reweight_name(base: str, label: str) -> str:
    """Return the expected cosmo_reweight fits filename for a label."""
    # base is e.g. 'IGWN-GWTC3p0-v3-GW191103_012549_PEDataRelease_mixed_cosmo.h5'
    stem = base[:-3]  # strip .h5
    return f"{stem}_reweight_{label}.fits"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-dir", required=True,
                        help="Directory containing *_cosmo.h5 files and skymaps")
    parser.add_argument("--ligo-skymap", required=True,
                        help="Path to ligo-skymap-from-samples executable")
    parser.add_argument("--jobs", type=int, default=16,
                        help="Number of parallel threads for ligo-skymap (-j)")
    parser.add_argument("--maxpts", type=int, default=2000,
                        help="Maximum posterior samples for skymap (--maxpts)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print commands without executing them")
    args = parser.parse_args()

    release_dir = os.path.abspath(args.release_dir)
    h5_files = sorted(
        f for f in os.listdir(release_dir) if f.endswith("_cosmo.h5")
    )

    if not h5_files:
        print(f"No *_cosmo.h5 files found in {release_dir}", file=sys.stderr)
        sys.exit(1)

    missing = []
    for h5_name in h5_files:
        h5_path = os.path.join(release_dir, h5_name)
        labels = get_h5_labels(h5_path)
        for label in labels:
            fits_name = cosmo_reweight_name(h5_name, label)
            fits_path = os.path.join(release_dir, fits_name)
            if not os.path.exists(fits_path):
                missing.append((h5_path, label, fits_path))

    if not missing:
        print("All cosmo_reweight skymaps already exist.")
        return

    print(f"Generating {len(missing)} missing skymap(s)...", file=sys.stderr)

    for h5_path, label, fits_path in missing:
        fits_name = os.path.basename(fits_path)
        cmd = [
            args.ligo_skymap,
            h5_path,
            "--path", f"{label}/posterior_samples",
            "--outdir", release_dir,
            "--fitsoutname", fits_name,
            "-j", str(args.jobs),
            "--maxpts", str(args.maxpts),
        ]
        print(" \\\n  ".join(cmd))
        if not args.dry_run:
            result = subprocess.run(cmd)
            if result.returncode != 0:
                print(f"ERROR: failed for {label} in {os.path.basename(h5_path)}",
                      file=sys.stderr)
                sys.exit(result.returncode)
        print()


if __name__ == "__main__":
    main()
