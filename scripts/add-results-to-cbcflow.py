"""
Add combined parameter estimation results to the CBCFlow library
"""

import os
import subprocess
import h5py as hdf5
import cbcflow
import sys
import yaml
from gpstime import gpstime
import click

def get_approved_analyses(event_name, ledger_path):
    """
    Get a list of approved bilby analyses for a given event from the ledger file.
    """
    with open(ledger_path, "r") as ledger_file:
        data = yaml.safe_load(ledger_file)
        
    events = {event['name']: event for event in data['events']}
    event = events[event_name]
    productions = {list(p.keys())[0]: list(p.values())[0] for p in event['productions']}
    
    approved_analyses = []
    for production, settings in productions.items():
        if (settings.get("review", None) is not None) and (len(settings['review'])>0 and settings['review'][-1]['status'].upper() == "APPROVED"):
            if settings['pipeline'] == "bilby":
                approximant = settings['waveform']['approximant']
                if approximant == "IMRPhenomXPHM":
                    approximant = "IMRPhenomXPHM-SpinTaylor"
                approved_analyses.append(approximant)
    return approved_analyses

@click.group()
def release():
    pass


@release.command("waveforms")
@click.option('--ledger', type=str, required=True, help='Path to the ledger.yml file')
@click.option('--version', type=str, default='v1', help='Version tag for the results to be checked')
@click.option('--event', type=str, required=True, help='Event name')
@click.option('--filename', type=str, default='C00-Combined.h5', help='Filename of the combined results file')
@click.option('--extras', type=str, help='Check for extra datasets not in the ledger')
def check_waveforms(event, ledger, version, filename, extras):    
    """
    Check that all available waveforms have been included in the combined results file for a given event.
    """
    ledger_waveforms = get_approved_analyses(event, ledger)
    print(ledger_waveforms)
    release_file = os.path.join("releases", version, event, filename)

    try:
        with hdf5.File(release_file) as datafile:
            waveforms = list(datafile.keys())
            waveforms.remove("history")
            waveforms.remove("version")

    except Exception as e:
        print(e)
        print(release_file, event)

    file_waveforms = []
    for waveform in waveforms:

        waveform = waveform.replace("C00:", "")
        waveform = waveform.replace(":", "-")

        file_waveforms.append(waveform)

    extras = extras.split(",") if extras is not None else []
    if (set(ledger_waveforms + extras)) - set(file_waveforms) != set():
        click.echo(f"Discrepancy for {event} in {version}")
        click.echo((set(ledger_waveforms + extras)) - set(file_waveforms))
        raise ValueError(f"Discrepancy for {event} in {version}")
    else:
        return True

@release.command("rename_waveforms")
@click.option('--filename', type=str, required=True, help='Path to the combined results file')
def rename_waveforms(filename):
    """
    Rename datasets in the combined results file to match datarelease naming conventions.
    """
    release_file = filename

    try:
        with hdf5.File(release_file) as datafile:
            waveforms = list(datafile.keys())
            waveforms.remove("history")
            waveforms.remove("version")
    except Exception as e:      
        print(e)

    for waveform in waveforms:

        new_name = waveform.

        subprocess.run(["summarymodify",
        "--delimiter=~",
        "--overwrite",
        "--labels", f"{waveform}~{new_name}",
        "--webdir", "tmp",
        "--samples", filename])

@release.command("rename_mixed")
@click.option('--filename', type=str, required=True, help='Path to the combined results file')
def rename_mixed(filename):
    """
    Rename datasets which have been combined with the summarycombine_posteriors command to "C00:Mixed"
    """
    release_file = filename

    try:
        with hdf5.File(release_file) as datafile:
            waveforms = list(datafile.keys())
            waveforms.remove("history")
            waveforms.remove("version")
    except Exception as e:      
        print(e)

    mixed = []
    for waveform in waveforms:
        if "_combined" in waveform:
            mixed.append(waveform)

    for waveform in mixed:
        if "XO4a" in waveform:
            new_name = "C00:Mixed+XO4a"
        else:   
            new_name = "C00:Mixed"

        subprocess.run(["summarymodify",
        "--delimiter=~",
        "--overwrite",
        "--labels", f"{waveform}~{new_name}",
        "--webdir", "tmp",
        "--samples", filename])

@release.command() 
@click.option('--ledger', type=str, required=True, help='Path to the ledger.yml file')
@click.option('--library', type=str, required=True, help='Path to the CBCFlow library')
@click.option('--version', type=str, default='v1', help='Version tag for the results to be added')
@click.option('--filename', type=str, default='C00-Combined.h5', help='Filename of the combined results file')
def add_skymaps():
    """
    Add skymaps to the CBCFlow library.
    """

    with open(ledger, "r") as datafile:
        data = yaml.safe_load(datafile)

    events = {event['name']: gpstime.fromgps(int(event['event time'])).strftime("GW%y%m%d_%H%M%S") for event in data['events']}

    library = cbcflow.core.database.LocalLibraryDatabase(library)
    
    for event, phone in events.items():

        release_file = os.path.join("releases", version, event, filename)

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
            analysis_output["UID"] = f"{version}-skymap-{waveform}"

            if analysis_output["UID"] in metadata_pe_results_uids:
                corresponding_analysis = metadata_pe_results[
                    metadata_pe_results_uids.index(analysis_output["UID"])
                ]
                continue
            else:
                corresponding_analysis = None

            analysis_output["InferenceSoftware"] = "PE Skymap Postprocessing"

            analysis_output["Notes"] = []



            analysis_output["Notes"].append(f"v3.0 Skymap for C00:{waveform}")
            analysis_output["PESummaryResultFile"] = {}
            analysis_output["PESummaryResultFile"]["Path"] = f"/home/pe.o4/GWTC4-fogg/releases/v4.0/{event}/C00-{waveform}.fits.gz"


            if os.path.exists(analysis_output["PESummaryResultFile"]["Path"]):

                pe.append(analysis_output)
                metadata.update(output)


                metadata.write_to_library(
                    message=f"Adding waveform to datarelease {version}", 
                    branch_name=f"add-{version}-skymap"
                )

            else:
                print(f"No skymap for {event} {waveform}")




@release.command("add_datafile") 
@click.option('--ledger', type=str, required=True, help='Path to the ledger.yml file')
@click.option('--library', type=str, required=True, help='Path to the CBCFlow library')
@click.option('--version', type=str, default='v1', help='Version tag for the results to be added')
@click.option('--filename', type=str, default='C00-Combined.h5', help='Filename of the combined results file')
def add_datafile(ledger, library, version, filename):
    """
    Add the combined results file to the CBCFlow library.
    """

    with open(ledger, "r") as datafile:
        data = yaml.safe_load(datafile)

    events = {event['name']: gpstime.fromgps(int(event['event time'])).strftime("GW%y%m%d_%H%M%S") for event in data['events']}

    library = cbcflow.core.database.LocalLibraryDatabase(library)

    for event, phone in events.items():

        release_file = os.path.join("releases", version, event, filename)

        if check_waveforms.callback(event, ledger, version, filename, extras="Mixed") is not True:
            print(f"Skipping {event} due to waveform discrepancy")
            continue


        try:
        
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
            analysis_output["UID"] = f"{version}-combined-results"

            if analysis_output["UID"] in metadata_pe_results_uids:
                corresponding_analysis = metadata_pe_results[
                    metadata_pe_results_uids.index(analysis_output["UID"])
                ]
            else:
                corresponding_analysis = None

            analysis_output["InferenceSoftware"] = "PE Postprocessing"

            analysis_output["Notes"] = []

            analysis_output["Notes"].append(f"{version} combined data release metafile")
            analysis_output["PESummaryResultFile"] = {}
            analysis_output["PESummaryResultFile"]["Path"] = os.path.abspath(release_file)

            pe.append(analysis_output)
            metadata.update(output)

            metadata.write_to_library(
                message=f"Adding PE datarelease {version}", 
                branch_name=f"add-{version}"
            )

        except Exception as e:
            print(e)
            print(event, phone)


if __name__ == "__main__":
    release()