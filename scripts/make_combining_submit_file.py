"""
Generate an HTCondor submit file for Stage 0 (review combined summary pages).

Each event gets its own job running:
    make ${WEB_PATH}/{event}/combined/samples/posterior_samples.h5

Run from the project directory:
    python ../scripts/make_combining_submit_file.py > ../condor/review_summarypages.sub

Pass an optional event name to generate a submit file for a single event:
    python ../scripts/make_combining_submit_file.py S230601bf

Set FORCE=1 in the environment to unconditionally rebuild all events (passes
-B to make and skips the "no samples on disk" guard):
    FORCE=1 python ../scripts/make_combining_submit_file.py > ../condor/review_summarypages.sub
"""
import sys
import os
sys.path.insert(0, "..")
from scripts.ledger import load_ledger, get_events, not_rejected_productions


def any_samples_exist(productions, web_path, event_name):
    """Return True if at least one pesummary sample file exists on disk."""
    return any(
        os.path.exists(
            f"{web_path}/{event_name}/{production}/pesummary/samples/posterior_samples.h5"
        )
        for production, _ in productions
    )

data = load_ledger(".asimov/ledger.yml")
events = get_events(data)

# The script is invoked after `cd $(PROJECT)` in the Makefile rule, so
# os.getcwd() is the project/ subdirectory.  The sub-makefiles and condor/
# directory live one level up (the Makefile root).
project_dir = os.path.dirname(os.getcwd())
web_path = os.environ.get("WEB_PATH", "/home/pe.o4/public_html/gwtc4-fogg")
env = os.environ.get("ENV", "/home/pe.o4/GWTC4-fogg/environment")
user = os.environ.get("ACCOUNTING_USER") or os.environ.get("USER", "pe.o4")
force = os.environ.get("FORCE", "").strip() in ("1", "true", "yes")
force_flag = "-B " if force else ""
ncpus = int(os.environ.get("NCPUS", 16))

description = """\
NCPUS = {ncpus}
executable = /usr/bin/make
initialdir = {project_dir}
arguments = {force_flag}-f combined_summary.mk SUMMARYPAGES={env}/bin/summarypages WEB_PATH={web_path} NCPUS={ncpus} {web_path}/{event}/combined/samples/posterior_samples.h5
request_disk = 20GB
request_memory = 16384MB
request_cpus = $(NCPUS)
getenv = CONDA_EXE,USER,LAL*,PATH,HOME
accounting_group_user = {user}
accounting_group = ligo.prod.o4.cbc.pe.bilby
log = condor/review_summarypages_{event}.log
output = condor/review_summarypages_{event}.out
error = condor/review_summarypages_{event}.err
batch_name = review-summarypages/{event}
queue

"""

filter_event = sys.argv[1] if len(sys.argv) > 1 else None

os.makedirs("condor", exist_ok=True)

for name, event in events.items():
    if filter_event and filter_event != name:
        continue
    productions = list(not_rejected_productions(event))
    if not productions:
        continue
    if not force and not any_samples_exist(productions, web_path, name):
        continue
    print(description.format(
        force_flag=force_flag,
        ncpus=ncpus,
        env=env,
        project_dir=project_dir,
        web_path=web_path,
        event=name,
        user=user,
    ))
