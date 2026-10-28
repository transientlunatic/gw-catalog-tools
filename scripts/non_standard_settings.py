#!/usr/bin/env python3
"""
Generate a markdown table of events whose settings differ from the defaults.

Reports events where:
- Any frame type differs from the global defaults (indicating deglitched frames)
- The minimum frequency differs from the standard 20 Hz for any active detector
"""

import argparse

from ledger import load_ledger

DEFAULT_MIN_FREQ = 20
ACTIVE_IFOS = ("H1", "L1", "V1")


def find_non_standard_events(ledger):
    """Yield (event_name, changed_frames, non_default_min_freqs) for non-standard events."""
    default_frames = ledger["data"]["frame types"]

    for event in ledger["events"]:
        name = event["name"]
        event_frames = event.get("data", {}).get("frame types", {})
        min_freqs = event.get("quality", {}).get("minimum frequency", {})

        changed_frames = {
            ifo: ft
            for ifo, ft in event_frames.items()
            if ft != default_frames.get(ifo)
        }

        changed_freqs = {
            ifo: f
            for ifo, f in min_freqs.items()
            if ifo in ACTIVE_IFOS and f != DEFAULT_MIN_FREQ
        }

        if changed_frames or changed_freqs:
            yield name, changed_frames, changed_freqs


def build_table(ledger):
    """Return a markdown table string."""
    rows = list(find_non_standard_events(ledger))
    if not rows:
        return "No events with non-standard settings found."

    lines = [
        "| Event | Detector | Deglitched frame type | f_low (Hz) |",
        "|-------|----------|-----------------------|------------|",
    ]
    for name, changed_frames, changed_freqs in rows:
        # Collect all detectors that have any non-standard setting
        all_ifos = sorted(set(changed_frames) | set(changed_freqs))
        for ifo in all_ifos:
            frame_str = f"`{changed_frames[ifo]}`" if ifo in changed_frames else "--"
            freq_str = str(changed_freqs[ifo]) if ifo in changed_freqs else "--"
            lines.append(f"| {name} | {ifo} | {frame_str} | {freq_str} |")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "ledger",
        nargs="?",
        default="ledger.yml",
        help="Path to the ledger YAML file (default: ledger.yml)",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Write output to a file instead of stdout",
    )
    args = parser.parse_args()

    ledger = load_ledger(args.ledger)
    table = build_table(ledger)

    if args.output:
        with open(args.output, "w") as f:
            f.write(table + "\n")
    else:
        print(table)


if __name__ == "__main__":
    main()
