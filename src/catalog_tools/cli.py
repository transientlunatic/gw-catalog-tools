"""The ``catalog-tools`` command line."""

import argparse
import json
import os
import sys

from . import __version__
from .batch import (DEFAULT_EXCLUDE, TEST_ORDER, discover_files, generate_markdown,
                    parse_path_ids, review_tree)
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


def write_markdown(data, path, arguments):
    """Write the markdown report to ``path``, or print it if there is none."""
    markdown = generate_markdown(data, emoji=not arguments.no_emoji, ignore=arguments.ignore)
    if path:
        with open(path, "w") as handle:
            handle.write(markdown + "\n")
        print(f"Markdown report written to {path}")
    else:
        print(markdown)


def batch(arguments):
    """Review many posterior samples files and write JSON and markdown summaries."""
    if arguments.h5files:
        files = [(*parse_path_ids(path), path) for path in arguments.h5files]
    elif arguments.root:
        try:
            files = list(discover_files(arguments.root, arguments.analysis, arguments.exclude))
        except FileNotFoundError as error:
            print(f"ERROR: {error}", file=sys.stderr)
            return 2
    else:
        print("ERROR: give H5 files or --root", file=sys.stderr)
        return 2
    if not files:
        print("No matching posterior_samples.h5 files found.", file=sys.stderr)
        return 1
    print(f"Reviewing {len(files)} analyses ...")

    data = review_tree(files, meta={"root": arguments.root,
                                    "analysis_patterns": arguments.analysis or [],
                                    "exclude": arguments.exclude})
    with open(arguments.output, "w") as handle:
        json.dump(data, handle, indent=2)
    print(f"Results saved to {arguments.output}")
    if arguments.markdown:
        write_markdown(data, arguments.markdown, arguments)

    records = [r for analyses in data["events"].values() for r in analyses.values()]
    return 1 if any(r["load_error"] or r["overall_status"] == "FAIL" for r in records) else 0


def report(arguments):
    """Render the JSON from ``batch`` as a markdown table."""
    with open(arguments.json_file) as handle:
        data = json.load(handle)
    write_markdown(data, arguments.output, arguments)
    return 0


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


def add_report_options(command):
    command.add_argument("--no-emoji", action="store_true", help="use PASS/WARN/FAIL rather than emoji")
    command.add_argument("--ignore", nargs="+", metavar="TEST",
                         help="tests to leave out of the markdown table (they stay in the JSON); "
                              f"one of {', '.join(TEST_ORDER)}")


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
        "batch", help="review many samples files and summarise them",
        description="Review every posterior_samples.h5 under a results root laid out as "
                    "ROOT/EVENT/ANALYSIS/posterior_samples.h5 (or the files named), write the "
                    "results as JSON and optionally a markdown table.  The exit status is 1 if "
                    "any analysis FAILs or can't be read, and 0 otherwise.")
    command.add_argument("h5files", nargs="*", metavar="H5FILE",
                         help="posterior_samples.h5 files (instead of --root)")
    command.add_argument("--root", "-r", metavar="DIR", help="results root; events are its subdirectories")
    command.add_argument("--analysis", "-a", nargs="+", metavar="PATTERN",
                         help="analysis names or globs to include, e.g. 'bilby-*' (default: all)")
    command.add_argument("--exclude", nargs="*", metavar="DIR", default=list(DEFAULT_EXCLUDE),
                         help=f"analysis directories to skip (default: {' '.join(DEFAULT_EXCLUDE)}); "
                              "give --exclude alone to skip none")
    command.add_argument("--output", "-o", metavar="JSON", default="review.json",
                         help="the JSON results file (default: review.json)")
    command.add_argument("--markdown", "-m", metavar="MD", help="also write a markdown report here")
    add_report_options(command)
    command.set_defaults(function=batch)

    command = commands.add_parser(
        "report", help="render the JSON from batch as a markdown table")
    command.add_argument("json_file", metavar="JSON", help="the file written by batch")
    command.add_argument("--output", "-o", metavar="MD", help="the markdown file (default: print it)")
    add_report_options(command)
    command.set_defaults(function=report)

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
