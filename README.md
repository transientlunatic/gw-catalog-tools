# gw-catalog-tools

Automated checks of gravitational-wave parameter estimation (PE) results, for
producing catalogues at scale.  Reviewing every posterior by eye does not scale
to a catalogue, so these checks flag the common failure modes — a posterior
railing against its prior, a signal which is really noise or a glitch, several
modes, a broken post-processing step — and leave the judgement to an analyst.

The checks work on PESummary `posterior_samples.h5` files.  They can be run from
the command line on one file, or as tests in an [asimov](https://asimov.readthedocs.io)
project, where they are run on each analysis awaiting review and write an HTML
report.

## Installation

```bash
pip install gw-catalog-tools            # the checks and the command line
pip install "gw-catalog-tools[asimov]"  # and what the asimov tests need
```

From a checkout, `pip install .` (or `pip install -e .`).  The checks need only
numpy, scipy and h5py.

## Command line

```bash
catalog-tools review posterior_samples.h5          # print a table of results
catalog-tools review posterior_samples.h5 --json   # machine-readable results
```

The exit status is 1 if any check FAILs, 2 if the file can't be read, and 0
otherwise.  Warnings don't affect it.

To review a whole tree of results laid out as
`ROOT/EVENT/ANALYSIS/posterior_samples.h5`, and summarise them in a table:

```bash
catalog-tools batch --root ROOT --analysis "bilby-*" -o review.json -m review.md
catalog-tools report review.json --no-emoji --ignore sky_localization
```

`batch` skips `generate-psd` directories unless `--exclude` says otherwise, and
exits 1 if any analysis FAILs or can't be read.  `report` renders the JSON again
without re-running anything.

## In an asimov project

```bash
cd my-project
catalog-tools install-tests tests     # writes tests/test_catalog_tools.py
python -m pytest tests
```

The module it writes is a single import of `catalog_tools.asimov_tests`, so
pytest finds the tests and upgrading the package upgrades them.  asimov must be
version 0.7 or later, with its bilby pipeline available.

By default the review tests run on uploaded bilby analyses which have not been
reviewed yet.  Environment variables change that:

| Variable | Effect |
|---|---|
| `PE_REVIEW_INCLUDE_REVIEWED=1` | also test analyses which are already approved, rejected or deprecated (useful for validating the checks against past decisions) |
| `PE_REVIEW_REPORT_DIR=DIR` | write the HTML reports to `DIR/<event>/<analysis>/automated_review.html` instead of beside the PESummary pages in the web root |

The samples file of an analysis is looked for in the project's web root
(`<webroot>/<event>/<analysis>/pesummary/samples/posterior_samples.h5`) and then
in its results store (`results/<event>/<analysis>/posterior_samples.h5`).

## The checks

A FAIL means the analysis should not be accepted as it is; a WARNING flags
something for an analyst to look at.  Thresholds are in `catalog_tools.review`.

| Check | Looks for | Result | Threshold |
|---|---|---|---|
| `chirp_mass_railing`, `distance_railing` | the prior cutting off real support: the end bin of a histogram of the posterior holds a significant fraction of the density of the highest bin | FAIL, with suggested new bounds | end bin > 2% of the highest of 50 bins (as in peconfigurator); not for the upper end of q, nor for the lower end of distance if the prior already extends below the suggested bound |
| `multimodality` | more than one significant mode in q, chirp mass or (log) distance, by counting peaks of a kernel density estimate | WARNING | prominence ≥ 8% of the highest peak and ≥ 3% of the probability |
| `diptest_unimodality` | Hartigan's dip test rejecting unimodality of q, chirp mass or distance; complements the peak count above | WARNING | p ≤ 0.01 |
| `snr_negative` | negative samples in a matched-filter or optimal SNR column | WARNING | any |
| `spin_prior_comparison` | spins pushed to the bounds, relative to the prior | WARNING | more than 10% more posterior than prior beyond 0.7 (a₁, a₂, χp) or ±0.9 (χeff) |
| `snr_plausibility` | a signal which isn't there, is too loud for a near-threshold candidate, or sits in one detector | WARNING | ln Bayes factor < 1 or network SNR < 4; network SNR > 20; one detector holding more than 75% of the SNR² |
| `effective_sample_size` | too few independent samples | FAIL / WARNING | the sampler's effective sample size < 1000 / < 2500 (the sample count if not stored) |
| `parameter_consistency` | broken post-processing: total mass, mass ratio and chirp mass inconsistent with the component masses | FAIL | relative error > 0.1% |
| `likelihood_sanity` | non-finite or constant log-likelihood or log-prior | FAIL | std of log-likelihood < 0.1 |
| `network_snr` | network SNR not the quadrature sum of the detector SNRs | FAIL | relative error > 1% |
| `sky_localization` | the 90% sky area (informational; large areas are normal near threshold) | INFO, or WARNING if collapsed to a couple of cells | |

The SNR and Bayes-factor thresholds are tuned for a catalogue of near-threshold
candidates, where a real signal has a network SNR of roughly 6–15.  Change the
constants in `catalog_tools.review` for other uses.

Two further tests of the asimov project's own setup are included: that the
sampler output hasn't lost samples (`TestSampleLoss`), and that the waveform
arguments in each analysis' configuration are the intended ones
(`TestWaveformConfiguration`).

### Limitations

* The multimodality check does not report a second mode with less than about
  10–20% of the probability, and a shoulder against a prior bound is not a mode;
  railing is the check which reports that.
* The checks screen; they do not replace looking at the data.  A glitch which
  produces a plausible posterior will pass them all.

## The rest of this repository

The repository also holds the scripts and Makefile which drive a catalogue
release from an asimov ledger.  They are not part of the installed package; see
[docs/release-pipeline.md](docs/release-pipeline.md).

## Development

```bash
pip install -e ".[test]"
pytest tests
```

The tests use synthetic data with known answers, and need neither asimov nor
a project.
