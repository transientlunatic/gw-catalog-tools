"""pytest configuration for catalog-tools release validation tests.

Add --ledger and --release-dir on the command line to point tests at a
specific Asimov project and release directory:

    pytest tests/ \\
        --ledger project/.asimov/ledger.yml \\
        --release-dir releases/v4.0

Any test that declares an ``event_name_and_data`` fixture argument is
automatically parametrized over every event in the ledger that has at
least one APPROVED bilby production, giving per-event test IDs in the
pytest output.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))


def pytest_addoption(parser):
    parser.addoption(
        "--ledger",
        default=".asimov/ledger.yml",
        help="Path to the Asimov project ledger.yml",
    )
    parser.addoption(
        "--release-dir",
        default="releases/staging",
        help="Path to the release output directory to validate",
    )


def pytest_generate_tests(metafunc):
    """Parametrize tests requesting event_name_and_data over approved ledger events."""
    if "event_name_and_data" not in metafunc.fixturenames:
        return
    ledger_path = metafunc.config.getoption("--ledger")
    try:
        from ledger import load_ledger, get_events, approved_bilby_productions
        data = load_ledger(ledger_path)
        events = get_events(data)
        params = [
            (name, event)
            for name, event in sorted(events.items())
            if list(approved_bilby_productions(event))
        ]
    except Exception:
        params = []
    metafunc.parametrize(
        "event_name_and_data",
        params,
        ids=[p[0] for p in params],
    )


@pytest.fixture(scope="session")
def release_dir(request):
    return request.config.getoption("--release-dir")


@pytest.fixture(scope="session")
def ledger_path(request):
    return request.config.getoption("--ledger")
