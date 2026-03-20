"""
Shared utilities for reading and querying the Asimov ledger.

All other scripts should import from here rather than re-implementing
ledger access or review-filter logic.
"""

import yaml
import numpy as np
import h5py


# ---------------------------------------------------------------------------
# Waveform label mapping
# ---------------------------------------------------------------------------

#: Map from approximant names stored in the ledger to canonical release labels.
RELEASE_LABELS = {
    "IMRPhenomXPHM": "IMRPhenomXPHM-SpinTaylor",
}


def release_label(approximant: str) -> str:
    """Return the canonical release label for an approximant name."""
    return RELEASE_LABELS.get(approximant, approximant)


# ---------------------------------------------------------------------------
# Ledger I/O
# ---------------------------------------------------------------------------

def load_ledger(path: str = ".asimov/ledger.yml") -> dict:
    """Load and return the raw ledger data."""
    with open(path, "r") as fh:
        return yaml.safe_load(fh)


def get_events(ledger: dict) -> dict:
    """Return a ``{name: event_dict}`` mapping from a loaded ledger."""
    return {event["name"]: event for event in ledger["events"]}


def get_productions(event: dict) -> dict:
    """Return a ``{production_name: settings_dict}`` mapping for one event."""
    return {
        list(p.keys())[0]: list(p.values())[0]
        for p in event["productions"]
    }


# ---------------------------------------------------------------------------
# Review filter predicates
# ---------------------------------------------------------------------------

def is_approved(settings: dict) -> bool:
    """
    Return True only if the production has an explicit APPROVED review.

    Unreviewed productions always return False.
    """
    review = settings.get("review") or []
    return len(review) > 0 and review[-1]["status"].upper() == "APPROVED"


def is_not_rejected(settings: dict) -> bool:
    """
    Return True if the production has no review, or its last review is not
    REJECTED.  Allows unreviewed productions through (used for Stage 0).
    """
    review = settings.get("review") or []
    if len(review) == 0:
        return True
    return review[-1]["status"].upper() != "REJECTED"


# ---------------------------------------------------------------------------
# Convenience iterators
# ---------------------------------------------------------------------------

def approved_bilby_productions(event: dict):
    """
    Yield ``(production_name, settings)`` for every APPROVED bilby production
    in *event*.
    """
    for name, settings in get_productions(event).items():
        if settings.get("pipeline") == "bilby" and is_approved(settings):
            yield name, settings


def not_rejected_bilby_productions(event: dict):
    """
    Yield ``(production_name, settings)`` for every non-rejected bilby
    production in *event* (includes unreviewed).
    """
    for name, settings in get_productions(event).items():
        if settings.get("pipeline") == "bilby" and is_not_rejected(settings):
            yield name, settings


# ---------------------------------------------------------------------------
# PESummary file helpers
# ---------------------------------------------------------------------------

def open_pesummary_samples(metafile: str):
    """
    Open a PESummary HDF5 file and return ``(posterior, prior)``.

    ``posterior`` is a structured numpy array; ``prior`` is a
    ``{param: array}`` dict.
    """
    with h5py.File(metafile, "r") as fh:
        analyses = [k for k in fh.keys() if k not in ("history", "version")]
        key = analyses[0]
        posterior = np.array(fh[key]["posterior_samples"])
        prior = {
            name: np.array(data)
            for name, data in fh[key]["priors"]["samples"].items()
        }
    return posterior, prior
