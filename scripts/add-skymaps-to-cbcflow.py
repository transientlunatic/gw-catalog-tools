import h5py as hdf5
import os
import cbcflow
import sys
import yaml
from gpstime import gpstime

with open("/home/pe.o4/GWTC4-fogg/project/.asimov/ledger.yml", "r") as datafile:
    data = yaml.safe_load(datafile)

events = {event['name']: gpstime.fromgps(int(event['event time'])).strftime("GW%y%m%d_%H%M%S") for event in data['events']}

print(events)

# events = ["S230601bf"]

library = cbcflow.core.database.LocalLibraryDatabase("/home/pe.o4/GWTC4/cbcflow-library")

# event = sys.argv[1]

for event, phone in events.items():
    release_file = f"/home/pe.o4/GWTC4-fogg/releases/v4.0/{event}/C00-Combined-XO4a.h5"

    try:
        with hdf5.File(release_file) as datafile:
            waveforms = list(datafile.keys())
            waveforms.remove("history")
            waveforms.remove("version")
    except Exception as e:
        print(e)
        print(release_file, event)
        continue

    for waveform in waveforms:

        waveform = waveform.replace("C00:", "")
        waveform = waveform.replace(":", "-")
        
        metadata = cbcflow.get_superevent(
            event, library=library
        )

        output = {}
        output["ParameterEstimation"] = {}
        pe = output["ParameterEstimation"]["Results"] = []


        metadata_pe_results = metadata["ParameterEstimation"]["Results"]

        metadata_pe_results_uids = cbcflow.core.utils.get_uids_from_object_array(
            metadata_pe_results
        )


        analysis_output = {}
        analysis_output["UID"] = f"v3-skymap-{waveform}"

        if analysis_output["UID"] in metadata_pe_results_uids:
            corresponding_analysis = metadata_pe_results[
                metadata_pe_results_uids.index(analysis_output["UID"])
            ]
            continue
        else:
            corresponding_analysis = None

        #analysis_output["InferenceSoftware"] = "PE Postprocessing"

        analysis_output["Notes"] = []



        analysis_output["Notes"].append(f"v3.0 Skymap for C00:{waveform}")
        analysis_output["PESummaryResultFile"] = {}
        analysis_output["PESummaryResultFile"]["Path"] = f"/home/pe.o4/GWTC4-fogg/releases/v4.0/{event}/C00-{waveform}.fits.gz"


        if os.path.exists(analysis_output["PESummaryResultFile"]["Path"]):

            pe.append(analysis_output)
            metadata.update(output)


            metadata.write_to_library(
                message="Adding waveform to datarelease v3", branch_name="add-reweighted"
            )

        else:
            print(f"No skymap for {event} {waveform}")
