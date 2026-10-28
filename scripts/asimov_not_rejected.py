"""
List non-rejected bilby analyses for an event.

Usage:
    python asimov_not_rejected.py <event_name> [approximant]

Prints the names of all non-rejected bilby productions for the given event,
optionally filtered to a specific waveform approximant.
"""
import sys
sys.path.insert(0, "..")
from scripts.ledger import load_ledger, get_events, not_rejected_productions

if len(sys.argv) < 2:
    print(f"Usage: {sys.argv[0]} <event_name> [approximant]", file=sys.stderr)
    sys.exit(1)

data = load_ledger(".asimov/ledger.yml")
events = get_events(data)

if sys.argv[1] not in events:
    print(f"Event '{sys.argv[1]}' not found in ledger.", file=sys.stderr)
    sys.exit(1)

event = events[sys.argv[1]]
approximant_filter = sys.argv[2] if len(sys.argv) > 2 else None

for production, settings in not_rejected_productions(event):
    if approximant_filter is None or settings["waveform"]["approximant"] == approximant_filter:
        print(production, end=" ")
print()
