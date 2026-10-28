"""
Stage 2A: generate Makefile rules to produce the release-quality combined
summary pages.

Unlike Stage 0, these pages:
  - Include only APPROVED bilby analyses
  - Apply full spin evolution (forwards and backwards)
  - Write to ${WEB_PATH}/{event}/combined-release/ so they never collide
    with the Stage 0 review pages

The combined metafile produced here is the input for Stage 2B (mixing).

Run from the project directory:
    python ../scripts/generate_release_summarypages.py > ../release_summary.mk
"""
import sys
sys.path.insert(0, "..")
from scripts.ledger import (
    load_ledger, get_events, approved_productions
)

data = load_ledger(".asimov/ledger.yml")
events = get_events(data)

# Emit .ONESHELL so that multi-line recipe loops work when this sub-makefile
# is invoked directly with -f (without the main Makefile's .ONESHELL in scope).
print(".ONESHELL:")
print("SHELL=/usr/bin/bash")
print()

all_targets = []

for name, event in events.items():
    productions = list(approved_productions(event))
    if not productions:
        continue

    labels = [p for p, _ in productions]
    samples = [
        f"${{WEB_PATH}}/{name}/{p}/pesummary/samples/posterior_samples.h5"
        for p, _ in productions
    ]
    delta_f = 1 / event["data"]["segment length"]

    target = f"${{WEB_PATH}}/{name}/combined-release/samples/posterior_samples.h5"
    all_targets.append(target)

    # Use $(wildcard ...) so Make does not fail when individual analysis files
    # have not been produced yet — they are created by external HTCondor jobs.
    samples_list = " ".join(samples)
    pairs_str = " \\\n\t  ".join(f"{l}:{s}" for l, s in zip(labels, samples))

    print(f"\n{target}: $(wildcard {samples_list})")
    # With .ONESHELL the entire recipe runs in one bash invocation.
    # We dynamically filter to whichever (label, sample) pairs exist on disk so
    # that --labels and -s always have the same count.
    print('\t_labels=""')
    print('\t_samples=""')
    print(f"\tfor _pair in {pairs_str}; do")
    print('\t  _label=$${_pair%%:*}')
    print('\t  _sample=$${_pair##*:}')
    print('\t  if [ -f "$$_sample" ]; then')
    print('\t    _labels="$$_labels $$_label"')
    print('\t    _samples="$$_samples $$_sample"')
    print('\t  fi')
    print('\tdone')
    print(f'\t[ -n "$$_labels" ] || (echo "No completed analyses for {name}"; exit 0)')
    print(f'\t${{SUMMARYPAGES}} -w $(dir $@)../ \\')
    print('\t--labels $$_labels \\')
    print('\t-s $$_samples \\')
    print('\t--multi_process $${NCPUS:-16} \\')
    print('\t--nsamples_for_skymap 2000 \\')
    print('\t--gw \\')
    print('\t--cosmology Planck15_lal \\')
    print('\t--redshift_method exact \\')
    print('\t--multi_threading_for_skymap 4 \\')
    print('\t--calculate_multipole_snr --calculate_precessing_snr \\')
    print('\t--evolve_spins \\')
    print('\t--evolve_spins_forwards \\')
    print('\t--evolve_spins_backwards hybrid_orbit_averaged \\')
    print(f'\t--delta_f {delta_f}')
    print()

if all_targets:
    print(".PHONY: all_release_pages\n")
    print("all_release_pages: " + " \\\n    ".join(all_targets) + "\n")
