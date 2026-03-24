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

recipe = """
${{WEB_PATH}}/{event}/combined-release/samples/posterior_samples.h5: {samples}"""

command = r"""	${{SUMMARYPAGES}} -w $(dir $@)../ \
--labels {labels} \
 -s $? \
 --multi_process 16 \
 --nsamples_for_skymap 2000 \
 --gw \
 --cosmology Planck15_lal \
 --redshift_method exact \
 --multi_threading_for_skymap 4 \
 --calculate_multipole_snr --calculate_precessing_snr \
 --evolve_spins \
 --evolve_spins_forwards \
 --evolve_spins_backwards hybrid_orbit_averaged \
 --delta_f {delta_f}
"""

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

    all_targets.append(
        f"${{WEB_PATH}}/{name}/combined-release/samples/posterior_samples.h5"
    )
    print(recipe.format(event=name, samples=" ".join(samples)))
    print(command.format(labels=" ".join(labels), delta_f=delta_f))

if all_targets:
    print(".PHONY: all_release_pages\n")
    print("all_release_pages: " + " \\\n    ".join(all_targets) + "\n")
