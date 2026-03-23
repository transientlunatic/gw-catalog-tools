"""
Stage 2B: generate Makefile rules to assemble the full combined release file.

For each event, takes the Stage 2A summarypages output (which has spin
evolution applied) and:
  1. Copies the combined metafile into releases/$(RELEASE)/{event}/C00-Combined.h5
  2. Renames each internal waveform label to the canonical C00:<WaveformName> form
  3. Runs summarycombine_posteriors to produce the C00:Mixed combined posterior

Only APPROVED bilby analyses are included.

Run from the project directory:
    python ../scripts/generate_mixing_command.py > ../mixing_release.mk
"""
import sys
sys.path.insert(0, "..")
from scripts.ledger import (
    load_ledger, get_events, approved_productions, release_label
)

data = load_ledger(".asimov/ledger.yml")
events = get_events(data)

# Stage 2A output lives alongside the Stage 0 review pages but in a separate
# subdirectory so the two can never contaminate each other.
recipe = """
releases/${{RELEASE}}/{event}/C00-Combined.h5: ${{WEB_PATH}}/{event}/combined-release/samples/posterior_samples.h5
\tmkdir -p $(@D)
\tcp $< $@"""

rename = """\t${{SUMMARYMODIFY}} \\
\t--samples $@ \\
\t--delimiter _ \\
\t--labels {production}_{new_label} \\
\t--overwrite \\
\t--webdir $(@D)
"""

combine = """\t${{SUMMARY_COMBINE_POSTERIORS}} \\
\t--samples $@ \\
\t--labels {labels} \\
\t--outdir $(@D) \\
\t--add_to_existing

\t${{SUMMARYMODIFY}} \\
\t--samples $(@D)/C00-Combined_combined.h5 \\
\t--delimiter / \\
\t--overwrite \\
\t--labels {labels_joined}_combined/C00:Mixed \\
\t--webdir $(@D)

\tmv $(@D)/C00-Combined_combined.h5 $@
"""

all_targets = []

for name, event in events.items():
    productions = list(approved_productions(event))
    if not productions:
        continue

    new_labels = []
    for production, settings in productions:
        approx = settings["waveform"]["approximant"]
        new_labels.append(f"C00:{release_label(approx)}")

    all_targets.append(f"releases/$(RELEASE)/{name}/C00-Combined.h5")
    print(recipe.format(event=name))
    for (production, _), new_label in zip(productions, new_labels):
        print(rename.format(production=production, new_label=new_label))

    if len(new_labels) > 1:
        print(combine.format(
            labels=" ".join(new_labels),
            labels_joined="_".join(new_labels),
        ))

if all_targets:
    print(".PHONY: all_combined_release\n")
    print("all_combined_release: " + " \\\n    ".join(all_targets) + "\n")
