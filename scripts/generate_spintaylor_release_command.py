"""
Stage 1: generate Makefile rules for the MVP (SpinTaylor-only) release.

For each event with exactly one APPROVED bilby IMRPhenomXPHM production,
a rule is emitted to copy the PESummary metafile into the release directory
and rename the internal label to 'C00:IMRPhenomXPHM-SpinTaylor'.

The source file (${WEB_PATH}/{event}/{production}/pesummary/samples/posterior_samples.h5)
must already be a PESummary metafile — i.e. the summarypages job must have
completed for that production.  Events where summarypages has not yet run
will fail at the summarymodify step; run summarypages first.

If an event has more than one qualifying production a warning is issued
to stderr and no rule is emitted for that event.  Unreviewed productions
are never included.

Run from the project directory:
    python ../scripts/generate_spintaylor_release_command.py > ../spintaylor_release.mk

Stage 1 files land under releases/$(RELEASE)/spintaylor/ to avoid any
conflict with the Stage 2B combined release at releases/$(RELEASE)/.
"""
import sys
sys.path.insert(0, "..")
from scripts.ledger import load_ledger, get_events, approved_productions

data = load_ledger(".asimov/ledger.yml")
events = get_events(data)

recipe = """
releases/${{RELEASE}}/spintaylor/{event}/C00-Combined.h5: ${{WEB_PATH}}/{event}/{production}/pesummary/samples/posterior_samples.h5
\tmkdir -p $(@D)
\tcp $< $@.tmp.h5
\t${{SUMMARYMODIFY}} --samples $@.tmp.h5 --delimiter _ --labels {production}_C00:IMRPhenomXPHM-SpinTaylor --overwrite --webdir /tmp
\tmv $@.tmp.h5 $@
"""

all_targets = []

for name, event in events.items():
    approved = [
        production
        for production, settings in approved_productions(event)
        if settings["waveform"]["approximant"] == "IMRPhenomXPHM"
    ]

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
    all_targets.append(f"releases/$(RELEASE)/spintaylor/{name}/C00-Combined.h5")
    print(recipe.format(event=name, production=production))

if all_targets:
    print(".PHONY: all_spintaylor\n")
    print("all_spintaylor: " + " \\\n    ".join(all_targets) + "\n")
