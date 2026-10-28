# catalog-tools release pipeline

> The automated review checks have moved into the installable `gw-catalog-tools` package (`src/catalog_tools`); see the top-level README. This page describes the repo-only release scripts.

Tools for assembling and releasing gravitational-wave parameter estimation (PE) data catalogs from [Asimov](https://asimov.docs.ligo.org)-managed projects.

## Overview

The workflow is driven by `$(PROJECT)/.asimov/ledger.yml` (default path `project/.asimov/ledger.yml`), which records every PE event, its individual analyses (productions), the pipeline and waveform used, and the review status of each analysis. Makefile rules are generated dynamically from the ledger so that the release always reflects the current state of the project.

All scripts that read the ledger import from the shared `scripts/ledger.py` module rather than re-implementing filter logic individually.

## Repository structure

```
Makefile                                   Top-level build orchestration
scripts/
  ledger.py                                Shared ledger I/O and review-filter utilities
  generate_combining_command.py            Stage 0: review summarypages make rules
  generate_release_summarypages.py         Stage 2A: release summarypages make rules
  generate_spintaylor_release_command.py   Stage 1: MVP SpinTaylor release make rules
  generate_mixing_command.py               Stage 2B: waveform mixing make rules
  skymaps.py                               Stage 2C: skymap make rules
  make_combining_submit_file.py            Condor submit file for Stage 0
  make_release_summarypages_submit.py      Condor submit file for Stage 2A
  make_mixing_submit_file.py               Condor submit file for Stage 2B
  gwtc4-comp.py                            Per-event waveform comparison quickview plot
  add-results-to-cbcflow.py               Add release files to the CBCFlow metadata library
  add-skymaps-to-cbcflow.py              Add skymaps to the CBCFlow metadata library
  event-table-to-blueprints.py            Convert a search summary table to Asimov blueprints
  asimov_not_rejected.py                  List non-rejected bilby analyses for an event
  conversions.py                          Chirp-mass / component-mass conversion utilities
  event-table-to-blueprints.py            Convert a search summary table to Asimov blueprints
  events.txt                              Full event list
  events-xo4a.txt                         Events with an IMRPhenomXO4a analysis
  events-nrsur.txt                        Events with an NRSur7dq4 analysis
  generate_release_report.py              Generate HTML + Markdown release summary report
tests/
  test_release_files.py                  Verify release file completeness and readability
  test_review.py                         Unit tests of the catalog_tools package (see the top-level README)
releases/                                  Release output directory (populated by make)
```

## Prerequisites

The following variables can be overridden on the command line or via the environment (defaults shown):

```bash
ENV             = /home/pe.o4/GWTC4-fogg/environment
WEB_PATH        = /home/pe.o4/public_html/gwtc4-fogg
RELEASE         = staging
PROJECT         = project          # path to the Asimov project directory
ACCOUNTING_USER =                  # HTCondor accounting user
BLUEPRINTS_DIR  = blueprints/analyses  # YAML blueprints applied by add-event
```

The following tools must be available inside `ENV`:

- `pesummary` (`summarypages`, `summarycombine_posteriors`, `summarymodify`)
- `ligo-skymap-from-samples`
- Python packages: `pyyaml`, `cbcflow`, `gpstime`, `h5py`, `click`, `numpy`

## Project bootstrap

If starting from scratch, the Makefile can initialise the Asimov project and
checkout the CBCFlow metadata library:

```bash
make project/.asimov/ledger.yml          # asimov init + apply defaults
make project/cbcflow-library             # clone CBCFlow library
```

The project init target applies the standard `production-pe`, `cbcflow`, and
`production-pe-priors` defaults from the [asimov data repository](https://git.ligo.org/asimov/data).
Override `PROJECT_NAME` or `ASIMOV_REGIME` if needed:

```bash
make project/.asimov/ledger.yml PROJECT_NAME="GWTC-5" ASIMOV_REGIME=testing
```

The CBCFlow clone requires two env vars — **never put credentials in the Makefile itself**:

```bash
export CBCFLOW_TOKEN=glpat-xxxx
export CBCFLOW_GIT_URL=git.ligo.org/cbc/projects/my-project.git
make project/cbcflow-library
```

## Event management

To add a single event and apply all blueprint YAML files from `BLUEPRINTS_DIR`:

```bash
make add-event-S230601bf
```

To add all events listed in `event_list.txt` (one event name per line):

```bash
make add-events
```

Blueprints in `BLUEPRINTS_DIR` are applied in shell glob order. Override
`BLUEPRINTS_DIR` to point at a different set of analysis recipes:

```bash
make add-event-S230601bf BLUEPRINTS_DIR=blueprints/rerun-analyses
```

## Release pipeline

The release is built in four stages. Stages 0 and 2A are compute-intensive and should be submitted via HTCondor. Stages 1, 2B, and 2C run quickly on the submit node.

### Stage 0 — Review combined summary pages

Runs `summarypages` across all non-rejected analyses per event using cheap settings (no spin evolution). Output lives at `${WEB_PATH}/{event}/combined/` and is used **only** for review — it is never part of any release. Regenerated as many times as needed during the review cycle.

```bash
make combined_summary.mk              # regenerate Stage 0 rules from ledger
make condor/review_summarypages.sub   # generate Condor submit file
condor_submit condor/review_summarypages.sub
```

### Stage 1 — MVP release (IMRPhenomXPHM-SpinTaylor only)

Copies the individual APPROVED pesummary file for each event into the release
directory and renames the internal label to `C00:IMRPhenomXPHM-SpinTaylor`.
The skymap is already present alongside the pesummary file and requires no
regeneration. Only analyses with an explicit **APPROVED** review status are
included; if an event has more than one qualifying analysis a warning is printed
and it is skipped.

Output lands under `releases/$(RELEASE)/spintaylor/` to keep it separate from
the Stage 2B combined release at `releases/$(RELEASE)/`. This avoids Makefile
target conflicts when both sub-makefiles are loaded simultaneously.

```bash
make spintaylor_release.mk            # regenerate Stage 1 rules from ledger
make all_spintaylor RELEASE=v4.0      # files → releases/v4.0/spintaylor/{event}/C00-Combined.h5
```

Once you have verified the output, submit to CBCFlow manually:

```bash
# Check waveforms are present (dry run):
python scripts/add-results-to-cbcflow.py waveforms \
    --ledger project/.asimov/ledger.yml \
    --release-dir releases/v4.0/spintaylor \
    --event S230601bf \
    --version v4.0

# Register files with CBCFlow:
python scripts/add-results-to-cbcflow.py add_datafile \
    --ledger project/.asimov/ledger.yml \
    --library project/cbcflow-library \
    --release-dir releases/v4.0/spintaylor \
    --version v4.0 \
    --spintaylor

# Register skymaps with CBCFlow:
python scripts/add-skymaps-to-cbcflow.py \
    --ledger project/.asimov/ledger.yml \
    --library project/cbcflow-library \
    --release-dir releases/v4.0/spintaylor \
    --version v4.0
```

### Stage 2A — Release-quality summarypages

Re-runs `summarypages` from the individual production pesummary files (APPROVED only) with full settings including spin evolution. Output goes to `${WEB_PATH}/{event}/combined-release/`, kept separate from Stage 0 to prevent unreviewed results contaminating the release. This is the input for Stage 2B.

```bash
make release_summary.mk                    # regenerate Stage 2A rules from ledger
make condor/release_summarypages.sub       # generate Condor submit file
condor_submit condor/release_summarypages.sub
```

### Stage 2B — Full combined release (all waveforms + Mixed)

Copies the Stage 2A combined metafile into `releases/$(RELEASE)/{event}/C00-Combined.h5`, renames each internal waveform label to the canonical `C00:<WaveformName>` form, then uses `summarycombine_posteriors` to create `C00:Mixed` combined posteriors.

```bash
make mixing_release.mk                  # regenerate Stage 2B rules from ledger
make condor/mixing_release.sub          # generate Condor submit file
condor_submit condor/mixing_release.sub
# or locally for one event:
make releases/v4.0/S230601bf/C00-Combined.h5 RELEASE=v4.0
```

### Stage 2C — Skymaps for full combined release

Generates a skymap for every waveform label in the combined HDF5, including `C00:Mixed`.

```bash
make skymaps.mk                    # write skymaps.mk from ledger
make all-skymaps RELEASE=v4.0
# or per-event:
make skymaps-S230601bf RELEASE=v4.0
```

`skymaps.mk` also contains rules for the Stage 1 SpinTaylor-only release
(one skymap per event, `C00-IMRPhenomXPHM-SpinTaylor.fits.gz`, under
`releases/$(RELEASE)/spintaylor/{event}/`):

```bash
make all-spintaylor-skymaps RELEASE=v4.0
# or per-event:
make spintaylor-skymap-S230601bf RELEASE=v4.0
```

Once you have verified the full combined release, submit to CBCFlow manually:

```bash
# Check waveforms are present (dry run):
python scripts/add-results-to-cbcflow.py waveforms \
    --ledger project/.asimov/ledger.yml \
    --event S230601bf \
    --version v4.0

# Register files with CBCFlow:
python scripts/add-results-to-cbcflow.py add_datafile \
    --ledger project/.asimov/ledger.yml \
    --library project/cbcflow-library \
    --version v4.0

# Register skymaps with CBCFlow:
python scripts/add-skymaps-to-cbcflow.py \
    --ledger project/.asimov/ledger.yml \
    --library project/cbcflow-library \
    --version v4.0
```

## Regenerating all sub-makefiles at once

The sub-makefiles (`combined_summary.mk`, `spintaylor_release.mk`, etc.) are
generated from the ledger and are only rebuilt automatically when the ledger is
**newer** than the existing `.mk` file. After pulling changes to the generator
scripts, always force a full regeneration before running a release build:

```bash
make makefiles
```

## Condor submit files

```bash
make condor-files
```

## Tests

```bash
make tests
# or directly:
cd project && python -m unittest discover tests
```

## Release verification

After the release HDF5 files and skymaps have been produced, verify that every
expected file exists, contains the correct waveform labels, and has at least one
posterior sample:

```bash
make release-check                    # check releases/staging/
make release-check RELEASE=v4.0      # or a specific release tag
```

This runs `tests/test_release_files.py`, which checks for each APPROVED event:

- `releases/{RELEASE}/{event}/C00-Combined.h5` exists and is non-empty
- Every expected `C00:{Waveform}` label (and `C00:Mixed` where applicable) is
  present in the HDF5
- Each label contains at least one posterior sample
- Each expected `C00-{Waveform}.fits.gz` skymap exists and is non-empty

If `releases/{RELEASE}/spintaylor/` exists it is checked automatically as well:
each event's `C00-Combined.h5` must exist and contain the
`C00:IMRPhenomXPHM-SpinTaylor` label with at least one sample.

Exits with code 0 on success and 1 on any failure.

## Release report

To generate a human-readable summary suitable for review or for pasting into a
Merge Request description:

```bash
make release-report                   # write to releases/staging/
make release-report RELEASE=v4.0
```

This runs `scripts/generate_release_report.py` and writes two files into the
release directory:

- `release_report.html` — a styled page with per-event waveform labels, sample
  counts, file sizes, and SHA-256 hashes
- `release_report.md` — the same information in GitHub-Flavoured Markdown,
  including a `sha256sums`-format block and per-event detail sections
  with skymap hashes, ready to be copied into an MR description

If `releases/{RELEASE}/spintaylor/` exists, both files automatically include a
**Stage 1 — SpinTaylor-only release** section covering that tree.

Both files are `.gitignore`d, as they are generated artefacts.
