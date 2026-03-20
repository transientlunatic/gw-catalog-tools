"""
This script takes an event name as an argument and prints a list of bilby analyses which are marked as "APPROVED" in the asimov project.
"""
import yaml
import sys

with open(".asimov/ledger.yml", "r") as ledger_file:
    data = yaml.safe_load(ledger_file)
    
events = {event['name']: event for event in data['events']}
event = events[sys.argv[1]]
productions = {list(p.keys())[0]: list(p.values())[0] for p in event['productions']}

for production, settings in productions.items():
    if (settings.get("review", None) is not None) or (len(settings['review'])>0 and settings['review'][-1]['status'].upper() != "REJECTED"):
        if len(sys.argv)==2:
            if settings['pipeline'] == "bilby":
                print(production, end=" ")
        elif settings['pipeline'] == "bilby" and settings['waveform']['approximant'] == sys.argv[2]:
            print(production, end=" ")
print()
