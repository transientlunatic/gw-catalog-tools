"""
Generate an HTCondor submit file for Stage 2B (mixing / combined release assembly).

Each event gets its own job running:
    make releases/{release}/{event}/C00-Combined.h5

The mixing step itself is fast, but it's submitted via Condor for consistency
and so it runs on a node with access to the shared filesystem.

Run from the project directory:
    python ../scripts/make_mixing_submit_file.py > ../condor/mixing_release.sub

Pass an optional event name to generate a submit file for a single event:
    python ../scripts/make_mixing_submit_file.py S230601bf
"""
import sys
import os
sys.path.insert(0, "..")
from scripts.ledger import load_ledger, get_events, approved_bilby_productions

data = load_ledger(".asimov/ledger.yml")
events = get_events(data)

makefile_path = os.environ.get("MAKEFILE_PATH", "/home/pe.o4/GWTC4-fogg/Makefile")
release = os.environ.get("RELEASE", "staging")
env = os.environ.get("ENV", "/home/pe.o4/GWTC4-fogg/environment")
user = os.environ.get("USER", "pe.o4")

description = """\
executable = {env}/bin/make
arguments = -f {makefile} releases/{release}/{event}/C00-Combined.h5
request_disk = 20GB
request_memory = 8192MB
request_cpus = 4
getenv = CONDA_EXE,USER,LAL*,PATH,HOME
accounting_group_user = {user}
accounting_group = ligo.prod.o4.cbc.pe.bilby
log = condor/mixing_{event}.log
output = condor/mixing_{event}.out
error = condor/mixing_{event}.err
batch_name = mixing/{event}
queue

"""

filter_event = sys.argv[1] if len(sys.argv) > 1 else None

os.makedirs("condor", exist_ok=True)

for name, event in events.items():
    if filter_event and filter_event != name:
        continue
    if not list(approved_bilby_productions(event)):
        continue
    print(description.format(
        env=env,
        makefile=makefile_path,
        release=release,
        event=name,
        user=user,
    ))
