"""
Generate Makefile rules for a data release containing only
APPROVED IMRPhenomXPHM-SpinTaylor samples.

For each event with exactly one APPROVED bilby IMRPhenomXPHM production,
a rule is emitted to copy the individual pesummary file into the release
directory and rename the internal label to 'C00:IMRPhenomXPHM-SpinTaylor'.

If an event has more than one qualifying production a warning is issued
to stderr and no rule is emitted for that event.  Unreviewed productions
are never included.
"""
import sys
import yaml

with open(".asimov/ledger.yml", "r") as ledger_file:
    data = yaml.safe_load(ledger_file)

events = {event['name']: event for event in data['events']}

recipe = """
releases/${{RELEASE}}/{event}/C00-Combined.h5: ${{WEB_PATH}}/{event}/{production}/pesummary/samples/posterior_samples.h5
\tmkdir -p $(@D)
\tcp $< $@"""

rename = """\t${{SUMMARYMODIFY}} \\
\t--samples $@ \\
\t--delimiter _ \\
\t--labels {production}_C00:IMRPhenomXPHM-SpinTaylor \\
\t--overwrite \\
\t--webdir $(@D)
"""

all_targets = []

for name, event in events.items():
    productions = {list(p.keys())[0]: list(p.values())[0] for p in event['productions']}

    approved = []
    for production, settings in productions.items():
        review = settings.get("review", None)
        if (review is not None
                and len(review) > 0
                and review[-1]['status'].upper() == "APPROVED"
                and settings['pipeline'] == "bilby"
                and settings['waveform']['approximant'] == "IMRPhenomXPHM"):
            approved.append(production)

    if len(approved) > 1:
        print(
            f"WARNING: {name} has {len(approved)} approved IMRPhenomXPHM productions "
            f"({', '.join(approved)}); skipping this event.",
            file=sys.stderr,
        )
        continue
    elif len(approved) == 0:
        continue

    production = approved[0]
    all_targets.append(f"releases/$(RELEASE)/{name}/C00-Combined.h5")
    print(recipe.format(event=name, production=production))
    print(rename.format(production=production))

if all_targets:
    print(".PHONY: all_spintaylor\n")
    print("all_spintaylor: " + " ".join(all_targets) + "\n")
