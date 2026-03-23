"""
Generate an HTCondor submit file for Stage 2A (release-quality summarypages).

Each event gets its own job running:
    make ${WEB_PATH}/{event}/combined-release/samples/posterior_samples.h5

These jobs are CPU- and memory-intensive because they apply spin evolution.

Run from the project directory:
    python ../scripts/make_release_summarypages_submit.py > ../condor/release_summarypages.sub

Pass an optional event name to generate a submit file for a single event:
    python ../scripts/make_release_summarypages_submit.py S230601bf
"""
import sys
import os
sys.path.insert(0, "..")
from scripts.ledger import load_ledger, get_events, approved_productions

data = load_ledger(".asimov/ledger.yml")
events = get_events(data)

makefile_path = os.environ.get("MAKEFILE_PATH", "/home/pe.o4/GWTC4-fogg/Makefile")
web_path = os.environ.get("WEB_PATH", "/home/pe.o4/public_html/gwtc4-fogg")
release = os.environ.get("RELEASE", "staging")

description = """\
executable = {env}/bin/make
arguments = -f {makefile} {web_path}/{event}/combined-release/samples/posterior_samples.h5
request_disk = 20GB
request_memory = 16384MB
request_cpus = 16
getenv = CONDA_EXE,USER,LAL*,PATH,HOME
accounting_group_user = {user}
accounting_group = ligo.prod.o4.cbc.pe.bilby
log = condor/release_summarypages_{event}.log
output = condor/release_summarypages_{event}.out
error = condor/release_summarypages_{event}.err
batch_name = release-summarypages/{event}
queue

"""

env = os.environ.get("ENV", "/home/pe.o4/GWTC4-fogg/environment")
user = os.environ.get("USER", "pe.o4")

filter_event = sys.argv[1] if len(sys.argv) > 1 else None

os.makedirs("condor", exist_ok=True)

for name, event in events.items():
    if filter_event and filter_event != name:
        continue
    if not list(approved_productions(event)):
        continue
    print(description.format(
        env=env,
        makefile=makefile_path,
        web_path=web_path,
        event=name,
        release=release,
        user=user,
    ))
