"""
This script takes an event name as an argument and prints a list of bilby analyses which are marked as "APPROVED" in the asimov project.
"""
import yaml
import sys

with open(".asimov/ledger.yml", "r") as ledger_file:
    data = yaml.safe_load(ledger_file)
    
events = {event['name']: event for event in data['events']}


description = """
executable = /home/pe.o4/GWTC4-fogg/environment/bin/make
arguments = -f /home/pe.o4/GWTC4-fogg/Makefile releases/staging/{event}/C00-Combined-XO4a.h5
request_disk = 20GB
request_memory = 8192MB
request_cpus = 16
getenv = CONDA_EXE,USER,LAL*,PATH,HOME
accounting_group_user = daniel.williams
accounting_group = ligo.prod.o4.cbc.pe.bilby
log = postprocessing-submit/xo4a_combined_{event}.log
output = postprocessing-submit/xo4a_combined_{event}.out
error = postprocessing-submit/xo4a_combined_{event}.err
batch_name = combined/{event}
queue
"""

if len(sys.argv) > 1:
    event_name = sys.argv[1]
else:
    event_name = None

for name, event in events.items():
    if event_name and event_name != name: continue
    # productions = {list(p.keys())[0]: list(p.values())[0] for p in event['productions']}
    # for production, settings in productions.items():
    #     if (settings.get("review", None) is not None) or (len(settings['review'])>0 and settings['review'][-1]['status'].upper() != "REJECTED"):
    #             if settings['pipeline'] == "bilby":
    #                 labels.append(production)
    #                 samples.append(f"${{WEB_PATH}}/{name}/{production}/pesummary/samples/posterior_samples.h5")

    print(description.format(event=name))
