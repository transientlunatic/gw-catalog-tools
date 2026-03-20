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
from scripts.ledger import load_ledger, get_events, not_rejected_bilby_productions

data = load_ledger(".asimov/ledger.yml")
events = get_events(data)

recipe = """
${{WEB_PATH}}/{event}/combined/samples/posterior_samples.h5: {samples}"""

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
 --delta_f {delta_f}
"""

all_targets = []

for name, event in events.items():
    labels = []
    samples = []
    for production, settings in not_rejected_bilby_productions(event):
        labels.append(production)
        samples.append(
            f"${{WEB_PATH}}/{name}/{production}/pesummary/samples/posterior_samples.h5"
        )

    if not labels:
        continue

    delta_f = 1 / event["likelihood"]["psd length"]
    all_targets.append(f"${{WEB_PATH}}/{name}/combined/samples/posterior_samples.h5")
    print(recipe.format(event=name, samples=" ".join(samples)))
    print(command.format(labels=" ".join(labels), delta_f=delta_f))

if all_targets:
    print(".PHONY: all_review_pages\n")
    print("all_review_pages: " + " \\\n    ".join(all_targets) + "\n")
