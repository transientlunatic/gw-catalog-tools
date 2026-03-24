# CLAUDE.md

Gravitational-wave parameter estimation (PE) catalog release tooling, built around [Asimov](https://asimov.docs.ligo.org)-managed projects.

## Build system

All orchestration goes through `make`. Sub-makefiles (`combined_summary.mk`, `release_summary.mk`, `spintaylor_release.mk`, `mixing_release.mk`, `skymaps.mk`) are **generated dynamically** from the Asimov ledger — never edit them by hand.

Regenerate all sub-makefiles at once:
```bash
make makefiles
```

Key configurable variables (override on the command line or via the environment):
- `ENV` — path to the Python environment containing pesummary tools
- `PROJECT` — path to the Asimov project directory (default: `project/`)
- `RELEASE` — release label (e.g. `v4.0`, default: `staging`)
- `WEB_PATH` — web output root
- `ACCOUNTING_USER` — HTCondor accounting user

## Shared ledger module

**All scripts must import from `scripts/ledger.py`** — never re-implement ledger access or review-filter logic inline. Key helpers:

- `load_ledger(path)` — load `ledger.yml`
- `get_events(ledger)` / `get_productions(event)` — dict accessors
- `is_approved(settings)` — True only for explicitly APPROVED productions
- `is_not_rejected(settings)` — True for unreviewed or non-REJECTED productions
- `approved_bilby_productions(event)` — iterator for Stage 1+

Waveform label canonicalisation also lives in `ledger.py` (`RELEASE_LABELS`, `release_label()`).

## Release pipeline stages

| Stage | Script | Condor? | Purpose |
|-------|--------|---------|---------|
| 0 | `generate_combining_command.py` | yes | Review summarypages (non-rejected, cheap settings) |
| 1 | `generate_spintaylor_release_command.py` | no | MVP SpinTaylor-only release (APPROVED only) |
| 2A | `generate_release_summarypages.py` | yes | Full-settings summarypages (APPROVED only) |
| 2B | `generate_mixing_command.py` | yes | Waveform mixing → `C00-Combined.h5` |
| 2C | `skymaps.py` | no | Skymaps for all labels in combined HDF5 |

Stage 0 output is **review-only** and must never contaminate a release.

## Tests

Tests live in `tests/` and use the `asimov` test framework:

```bash
# Run all tests
pytest tests/

# Key test files
tests/test_release_files.py       # release file completeness / readability
tests/test_automated_review.py    # railing, ESS, SNR consistency, etc.
tests/test_samples.py
tests/test_railing.py
```

## Credentials

CBCFlow operations need environment variables — never put credentials in the Makefile:
```bash
export CBCFLOW_TOKEN=glpat-xxxx
export CBCFLOW_GIT_URL=git.ligo.org/cbc/projects/my-project.git
```
