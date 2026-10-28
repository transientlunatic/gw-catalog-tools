"""
Stage 0: generate Makefile rules for combined review summary pages.

For each event, produces a summarypages job across all non-rejected bilby
analyses using cheap settings (no spin evolution) suitable for review.
Output lands in ${WEB_PATH}/{event}/combined/ and is NOT used in any release.

Run from the project directory:
    python ../scripts/generate_combining_command.py > ../combined_summary.mk
"""
import sys
sys.path.insert(0, "..")
from scripts.ledger import load_ledger, get_events, not_rejected_productions

data = load_ledger(".asimov/ledger.yml")
events = get_events(data)

# Emit .ONESHELL so that multi-line recipe loops work when this sub-makefile
# is invoked directly with -f (without the main Makefile's .ONESHELL in scope).
print(".ONESHELL:")
print("SHELL=/usr/bin/bash")
print()

all_targets = []

for name, event in events.items():
    labels = []
    samples = []
    for production, settings in not_rejected_productions(event):
        labels.append(production)
        samples.append(
            f"${{WEB_PATH}}/{name}/{production}/pesummary/samples/posterior_samples.h5"
        )

    if not labels:
        continue

    delta_f = 1 / event["data"]["segment length"]
    target = f"${{WEB_PATH}}/{name}/combined/samples/posterior_samples.h5"
    all_targets.append(target)

    # Use $(wildcard ...) so Make does not fail when individual analysis files
    # have not been produced yet — they are created by external HTCondor jobs.
    samples_list = " ".join(samples)
    # Shell loop pairs: "label:sample_path" entries separated by whitespace.
    # Each pair is on its own continuation line so very long events stay readable.
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
    print(f'\t--delta_f {delta_f}')
    print()

if all_targets:
    print(".PHONY: all_review_pages\n")
    print("all_review_pages: " + " \\\n    ".join(all_targets) + "\n")
