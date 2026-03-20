"""
This script takes an event name as an argument and prints a list of bilby analyses which are marked as "APPROVED" in the asimov project.
"""
import yaml
import sys

with open(".asimov/ledger.yml", "r") as ledger_file:
    data = yaml.safe_load(ledger_file)
    
events = {event['name']: event for event in data['events']}


recipe = """
releases/${{RELEASE}}/{event}/C00-Combined-XO4a.h5: ${{WEB_PATH}}/{event}/combined/samples/posterior_samples.h5
	mkdir -p $(@D)
	cp $< $@
"""

rename = r"""	${{SUMMARYMODIFY}} \
	--samples $@ \
	--delimiter _ \
	--labels {analysis}_{new_name} \
	--overwrite \
	--webdir $(@D)
"""

combine=r"""
	${{SUMMARY_COMBINE_POSTERIORS}} \
	--samples $@ \
	--labels {labels} \
	--outdir $(@D) \
	--add_to_existing

	${{SUMMARYMODIFY}} \
	--samples $(@D)/C00-Combined-XO4a_combined.h5 \
	--delimiter / \
	--overwrite \
	--labels {labels_join}_combined/C00:Mixed \
	--webdir $(@D)

	mv $(@D)/C00-Combined-XO4a_combined.h5 $@ 
"""

combine_xo4a=r"""
	${{SUMMARY_COMBINE_POSTERIORS}} \
	--samples $@ \
	--labels {labels} \
	--outdir $(@D) \
	--add_to_existing

	${{SUMMARYMODIFY}} \
	--samples $(@D)/C00-Combined-XO4a_combined.h5 \
	--delimiter / \
	--overwrite \
	--labels {labels_join}_combined/C00:Mixed+XO4a \
	--webdir $(@D)

	mv $(@D)/C00-Combined-XO4a_combined.h5 $@ 
"""

for name, event in events.items():
    labels = []
    waveforms = []
    productions = {list(p.keys())[0]: list(p.values())[0] for p in event['productions']}
    for production, settings in productions.items():
        if (settings.get("review", None) is None) or (settings.get("review", None) is {}) or (len(settings['review'])>0 and settings['review'][-1]['status'].upper() != "REJECTED"):
            if settings['pipeline'] == "bilby":
                labels.append(production)
                waveforms.append(settings['waveform']['approximant'])
                    
    print(recipe.format(event=name))
    new_names = []
    for label, waveform in zip(labels, waveforms):
        if waveform == "IMRPhenomXPHM":
            new_name = "C00:IMRPhenomXPHM-SpinTaylor"
        else:
            new_name = f"C00:{waveform}"
        new_names.append(new_name)
        print(rename.format(event=name, analysis=label, new_name=new_name))
    if "IMRPhenomXO4a" in waveforms:
        new_names_no_xo4a = new_names
        new_names_no_xo4a.remove("C00:IMRPhenomXO4a")
        print(combine.format(labels=" ".join(new_names), labels_join="_".join(new_names_no_xo4a)))
    if "IMRPhenomXO4a" in waveforms:
        print(combine_xo4a.format(labels=" ".join(new_names), labels_join="_".join(new_names)))
