"""
Add release skymaps to the CBCFlow metadata library.

For each event, registers every skymap FITS file present in the release
directory with CBCFlow.  Skips a waveform if the file does not exist on disk.

Usage
-----
    python add-skymaps-to-cbcflow.py \\
        --ledger project/.asimov/ledger.yml \\
        --library /path/to/cbcflow-library \\
        --version v4.0 \\
        --release-dir releases/v4.0
"""

import os
import sys
import click
import cbcflow

sys.path.insert(0, os.path.dirname(__file__))
from ledger import load_ledger, get_events, approved_productions, release_label


@click.command()
@click.option("--ledger", required=True, help="Path to the ledger.yml file")
@click.option("--library", required=True, help="Path to the CBCFlow library")
@click.option("--version", default="v1", show_default=True, help="Release version tag")
@click.option("--release-dir", default=None,
              help="Root directory of release files (default: releases/{version})")
def add_skymaps(ledger, library, version, release_dir):
    """Register skymap FITS files for all approved waveforms with CBCFlow."""
    from gpstime import gpstime as _gpstime
    import yaml

    if release_dir is None:
        release_dir = os.path.join("releases", version)

    with open(ledger) as fh:
        raw = yaml.safe_load(fh)

    events = {
        event["name"]: _gpstime.fromgps(int(event["event time"])).strftime("GW%y%m%d_%H%M%S")
        for event in raw["events"]
    }

    all_events = get_events(load_ledger(ledger))
    db = cbcflow.core.database.LocalLibraryDatabase(library)

    for event_name, gw_name in events.items():
        event = all_events[event_name]
        waveform_labels = [
            release_label(settings["waveform"]["approximant"])
            for _, settings in approved_productions(event)
        ]
        # Include the Mixed posterior skymap if multiple waveforms were combined
        if len(waveform_labels) > 1:
            waveform_labels.append("Mixed")

        for waveform in waveform_labels:
            fits_name = f"C00-{waveform.replace(':', '-')}.fits.gz"
            fits_path = os.path.join(release_dir, event_name, fits_name)

            if not os.path.exists(fits_path):
                click.echo(f"No skymap for {event_name} {waveform} — skipping.")
                continue

            uid = f"{version}-skymap-{waveform}"

            try:
                metadata = cbcflow.get_superevent(event_name, library=db)
                pe_results = metadata["ParameterEstimation"]["Results"]
                existing_uids = cbcflow.core.utils.get_uids_from_object_array(pe_results)

                if uid in existing_uids:
                    click.echo(f"{event_name} {uid} already in CBCFlow — skipping.")
                    continue

                output = {"ParameterEstimation": {"Results": []}}
                entry = {
                    "UID": uid,
                    "InferenceSoftware": "PE Skymap Postprocessing",
                    "Notes": [f"{version} skymap for C00:{waveform}"],
                    "PESummaryResultFile": {"Path": os.path.abspath(fits_path)},
                }
                output["ParameterEstimation"]["Results"].append(entry)
                metadata.update(output)
                metadata.write_to_library(
                    message=f"Adding {version} skymap for {waveform}",
                    branch_name=f"add-{version}-skymap",
                )
                click.echo(f"Added {event_name} {uid} to CBCFlow.")

            except Exception as exc:
                click.echo(f"Error processing {event_name} {waveform}: {exc}")


if __name__ == "__main__":
    add_skymaps()
