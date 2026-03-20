"""
Add parameter estimation results to the CBCFlow metadata library.

Commands
--------
waveforms    Check that all approved waveforms are present in a release file.
add_datafile Register a combined release HDF5 with CBCFlow.

Usage examples
--------------
Check waveforms:
    python add-results-to-cbcflow.py waveforms \\
        --ledger project/.asimov/ledger.yml \\
        --event S230601bf --version v4.0

Add a release file:
    python add-results-to-cbcflow.py add_datafile \\
        --ledger project/.asimov/ledger.yml \\
        --library /path/to/cbcflow-library \\
        --version v4.0
"""

import os
import subprocess
import sys
import h5py as hdf5
import cbcflow
import click

sys.path.insert(0, os.path.dirname(__file__))
from ledger import load_ledger, get_events, approved_bilby_productions, release_label


def _waveform_labels_in_file(release_file):
    """Return the list of waveform label strings from a release HDF5."""
    with hdf5.File(release_file) as fh:
        keys = [k for k in fh.keys() if k not in ("history", "version")]
    return [k.replace("C00:", "").replace(":", "-") for k in keys]


def _get_approved_waveform_labels(event_name, ledger_path):
    """Return canonical release labels for all approved waveforms of an event."""
    data = load_ledger(ledger_path)
    events = get_events(data)
    return [
        release_label(settings["waveform"]["approximant"])
        for _, settings in approved_bilby_productions(events[event_name])
    ]


@click.group()
def release():
    pass


@release.command("waveforms")
@click.option("--ledger", required=True, help="Path to the ledger.yml file")
@click.option("--version", default="v1", show_default=True, help="Release version tag")
@click.option("--event", required=True, help="Event name (e.g. S230601bf)")
@click.option("--filename", default="C00-Combined.h5", show_default=True,
              help="Release HDF5 filename")
@click.option("--extras", default=None,
              help="Comma-separated extra labels expected in the file (e.g. Mixed)")
def check_waveforms(event, ledger, version, filename, extras):
    """Check that all approved waveforms are present in the release file."""
    ledger_waveforms = _get_approved_waveform_labels(event, ledger)
    release_file = os.path.join("releases", version, event, filename)

    try:
        file_waveforms = _waveform_labels_in_file(release_file)
    except Exception as exc:
        raise click.ClickException(f"Could not read {release_file}: {exc}")

    extra_list = extras.split(",") if extras else []
    expected = set(ledger_waveforms + extra_list)
    missing = expected - set(file_waveforms)
    if missing:
        raise click.ClickException(
            f"Waveform discrepancy for {event} in {version}: missing {missing}"
        )
    click.echo(f"{event}: all expected waveforms present.")
    return True


@release.command("rename_mixed")
@click.option("--filename", required=True, help="Path to the combined release HDF5")
def rename_mixed(filename):
    """
    Rename combined-posterior datasets to C00:Mixed.

    Datasets whose label contains '_combined' (as produced by
    summarycombine_posteriors) are renamed to 'C00:Mixed'.
    """
    try:
        with hdf5.File(filename) as fh:
            waveforms = [k for k in fh.keys() if k not in ("history", "version")]
    except Exception as exc:
        raise click.ClickException(f"Could not read {filename}: {exc}")

    for waveform in waveforms:
        if "_combined" not in waveform:
            continue
        new_name = "C00:Mixed"
        subprocess.run(
            ["summarymodify",
             "--delimiter=~",
             "--overwrite",
             "--labels", f"{waveform}~{new_name}",
             "--webdir", "tmp",
             "--samples", filename],
            check=True,
        )


@release.command("add_datafile")
@click.option("--ledger", required=True, help="Path to the ledger.yml file")
@click.option("--library", required=True, help="Path to the CBCFlow library")
@click.option("--version", default="v1", show_default=True, help="Release version tag")
@click.option("--filename", default="C00-Combined.h5", show_default=True,
              help="Release HDF5 filename")
def add_datafile(ledger, library, version, filename):
    """Register the combined release HDF5 with CBCFlow for all events."""
    from gpstime import gpstime as _gpstime
    import yaml

    with open(ledger) as fh:
        data = yaml.safe_load(fh)

    events = {
        event["name"]: _gpstime.fromgps(int(event["event time"])).strftime("GW%y%m%d_%H%M%S")
        for event in data["events"]
    }

    db = cbcflow.core.database.LocalLibraryDatabase(library)

    for event_name, gw_name in events.items():
        release_file = os.path.join("releases", version, event_name, filename)

        # Validate before touching CBCFlow
        try:
            check_waveforms.callback(
                event=event_name, ledger=ledger, version=version,
                filename=filename, extras="Mixed",
            )
        except SystemExit:
            click.echo(f"Skipping {event_name}: waveform validation failed.")
            continue

        try:
            metadata = cbcflow.get_superevent(event_name, library=db)

            pe_results = metadata["ParameterEstimation"]["Results"]
            existing_uids = cbcflow.core.utils.get_uids_from_object_array(pe_results)
            uid = f"{version}-combined-results"

            output = {"ParameterEstimation": {"Results": []}}
            entry = {
                "UID": uid,
                "InferenceSoftware": "PE Postprocessing",
                "Notes": [f"{version} combined data release metafile"],
                "PESummaryResultFile": {"Path": os.path.abspath(release_file)},
            }
            output["ParameterEstimation"]["Results"].append(entry)
            metadata.update(output)
            metadata.write_to_library(
                message=f"Adding PE datarelease {version}",
                branch_name=f"add-{version}",
            )
            click.echo(f"Added {event_name} ({uid}) to CBCFlow.")

        except Exception as exc:
            click.echo(f"Error processing {event_name}: {exc}")


if __name__ == "__main__":
    release()
