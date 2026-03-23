"""
Stage 2C: generate Makefile rules to produce skymaps for the full combined
release.

A skymap is generated for every waveform label present in the combined
release HDF5 (including C00:Mixed). Rules are written to skymaps.mk.

Run from the project directory (after the release HDF5 files exist):
    python ../scripts/skymaps.py .asimov/ledger.yml
"""
import sys
sys.path.insert(0, "..")
from scripts.ledger import (
    load_ledger, get_events, approved_productions, release_label
)

if len(sys.argv) < 2:
    print(f"Usage: {sys.argv[0]} <ledger.yml>", file=sys.stderr)
    sys.exit(1)

data = load_ledger(sys.argv[1])
events = get_events(data)

all_skymap_targets = []

with open("skymaps.mk", "w") as skymap_file:

    for name, event in events.items():
        productions = list(approved_productions(event))
        if not productions:
            continue

        # The waveform labels that will be in the release HDF5, plus Mixed
        waveform_labels = [
            release_label(settings["waveform"]["approximant"])
            for _, settings in productions
        ]
        if len(waveform_labels) > 1:
            waveform_labels.append("Mixed")

        release_file = f"releases/$(RELEASE)/{name}/C00-Combined.h5"
        event_skymap_targets = []

        for waveform in waveform_labels:
            fits_name = f"C00-{waveform.replace(':', '-')}.fits.gz"
            target = f"releases/$(RELEASE)/{name}/{fits_name}"
            event_skymap_targets.append(target)
            all_skymap_targets.append(target)

            skymap_file.write(f"{target}: {release_file}\n")
            skymap_file.write(f"\tmkdir -p $(@D)\n")
            skymap_file.write(
                f"\tligo-skymap-from-samples $^ -j 32 --maxpts 2000 "
                f"--outdir $(@D) "
                f"--path 'C00:{waveform}/posterior_samples' "
                f"--fitsoutname {fits_name}\n"
            )
            skymap_file.write("\n")

        # Convenience phony target per event
        skymap_file.write(f".PHONY: skymaps-{name}\n")
        skymap_file.write(
            f"skymaps-{name}: {' '.join(event_skymap_targets)}\n\n"
        )

    if all_skymap_targets:
        skymap_file.write(".PHONY: all-skymaps\n")
        skymap_file.write(
            "all-skymaps: " + " \\\n    ".join(all_skymap_targets) + "\n"
        )

