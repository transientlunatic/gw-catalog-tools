import yaml
from gpstime import gpstime
import h5py as hdf5

import sys

with open(sys.argv[1], "r") as datafile:
    data = yaml.safe_load(datafile)

events = {event['name']: event for event in data['events']}

def waveforms_for_event(event):

    productions = {list(p.keys())[0]: list(p.values())[0] for p in events[event]['productions'] }
    out = []
    for production, settings in productions.items():
        if len(settings['review'])>0:
            if settings['review'][-1]['status'].upper() == "APPROVED" and settings['pipeline'] == "bilby":
                out.append(settings['waveform']['approximant'])
    return out


with open("skymaps.mk", "w") as skymap_file:

    events_in = []

    events["S230518h"] = "GW230518_125908"
    events["S230529ay"] = "GW230529_181500"
    
    for name, event in events.items():

        release_file = f"releases/staging/{name}/C00-Combined-XO4a.h5"

        try:
        
           with hdf5.File(release_file) as datafile:
               waveforms = list(datafile.keys())
               waveforms.remove("history")
               waveforms.remove("version")

        except FileNotFoundError:
            print(name)
            continue
               
        events_in += [name]
        for waveform in waveforms:
            waveform = waveform.replace("C00:", "")
            skymap_file.write(f"releases/${{RELEASE}}/{name}/C00-{waveform.replace(':', '-')}.fits.gz: {release_file} \n")
            skymap_file.write(f"\tmkdir -p releases/${{RELEASE}}/{name}/\n")
            skymap_file.write(f"\tligo-skymap-from-samples $^ -j 32 --maxpts 2000 --outdir $(@D) --path 'C00:{waveform}/posterior_samples' --fitsoutname C00-{waveform.replace(':','-')}.fits.gz\n")
            skymap_file.write("\n")

        skymap_file.write(f"skymaps/{name}: ")

        for waveform in waveforms:
            waveform = waveform.replace("C00:", "")
            skymap_file.write(f"releases/${{RELEASE}}/{name}/C00-{waveform.replace(':', '-')}.fits.gz ")
        skymap_file.write("\n")

    skymap_file.write(f"all-skymaps: skymaps/{' skymaps/'.join(events_in)}\n")
    skymap_file.write(".PHONY: all-skymaps\n")
