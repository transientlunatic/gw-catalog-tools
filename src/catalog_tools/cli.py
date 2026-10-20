"""The ``catalog-tools`` command line."""

import argparse
import json
import os
import sys

from . import __version__
from .review import open_pesummary_samples, read_run_info, run_all_checks

SHIM = '''"""Run the catalog-tools tests in this asimov project."""

from catalog_tools.asimov_tests import *  # noqa: F401,F403
'''


def review(arguments):
    """Run every check on one posterior samples file."""
    try:
        posterior, prior = open_pesummary_samples(arguments.h5file)
        info = read_run_info(arguments.h5file)
    except Exception as error:
        print(f"ERROR: failed to read {arguments.h5file}: {error}", file=sys.stderr)
        return 2

    results = run_all_checks(posterior, prior, info)

    if arguments.json:
        print(json.dumps(results, indent=2, default=float))
    else:
        print(f"{arguments.h5file}: {len(posterior)} samples, {len(posterior.dtype.names)} parameters\n")
        counts = {}
        for name, result in results.items():
            counts[result["status"]] = counts.get(result["status"], 0) + 1
            print(f"{result['status']:8} {name.replace('_', ' '):24} {result['message']}")
        print("\n" + ", ".join(f"{count} {status}" for status, count in sorted(counts.items())))

    return 1 if any(result["status"] == "FAIL" for result in results.values()) else 0


def install_tests(arguments):
    """Write a test module into an asimov project's tests directory."""
    os.makedirs(arguments.directory, exist_ok=True)
    path = os.path.join(arguments.directory, "test_catalog_tools.py")
    if os.path.exists(path) and not arguments.force:
        print(f"{path} already exists; use --force to replace it", file=sys.stderr)
        return 1
    with open(path, "w") as handle:
        handle.write(SHIM)
    print(f"Wrote {path}\nRun the tests from the project directory with:\n"
          f"    python -m pytest {arguments.directory}")
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="catalog-tools",
        description="Automated checks of gravitational-wave parameter estimation results.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)

    command = commands.add_parser(
        "review", help="run the checks on a posterior samples file",
        description="Run the checks on a PESummary posterior samples file.  The exit "
                    "status is 1 if any check FAILs, 2 if the file can't be read, and "
                    "0 otherwise; warnings do not affect it.")
    command.add_argument("h5file", help="path to posterior_samples.h5")
    command.add_argument("--json", action="store_true", help="print the results as JSON")
    command.set_defaults(function=review)

    command = commands.add_parser(
        "install-tests", help="add the asimov tests to a project",
        description="Write a test module into an asimov project's tests directory, "
                    "which imports the catalog-tools tests so that pytest finds them.")
    command.add_argument("directory", nargs="?", default="tests",
                         help="the project's tests directory (default: tests)")
    command.add_argument("--force", action="store_true", help="replace an existing module")
    command.set_defaults(function=install_tests)

    arguments = parser.parse_args(argv)
    return arguments.function(arguments)


if __name__ == "__main__":
    sys.exit(main())
