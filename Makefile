
.ONESHELL:
SHELL=/usr/bin/bash

# ---------------------------------------------------------------------------
# Configurable variables — override on the command line or in the environment
# ---------------------------------------------------------------------------
ENV             ?= /home/pe.o4/GWTC4-fogg/environment
WEB_PATH        ?= /home/pe.o4/public_html/gwtc4-fogg
RELEASE         ?= staging
PROJECT         ?= project
ACCOUNTING_USER ?=
# Directory that contains per-analysis blueprint YAML files applied by add-event
BLUEPRINTS_DIR  ?= blueprints/analyses

MAKEFILE_PATH ?= $(abspath $(lastword $(MAKEFILE_LIST)))

ASIMOV                    = ${ENV}/bin/asimov
SUMMARYPAGES              = ${ENV}/bin/summarypages
SUMMARY_COMBINE_POSTERIORS= ${ENV}/bin/summarycombine_posteriors
SUMMARYMODIFY             = ${ENV}/bin/summarymodify
LIGO_SKYMAP               = ${ENV}/bin/ligo-skymap-from-samples
PYTHON                    = ${ENV}/bin/python

export ENV WEB_PATH RELEASE PROJECT MAKEFILE_PATH ACCOUNTING_USER

# ---------------------------------------------------------------------------
# Generated sub-makefiles
# ---------------------------------------------------------------------------
include combined_summary.mk
include release_summary.mk
include spintaylor_release.mk
include mixing_release.mk
-include skymaps.mk

# ---------------------------------------------------------------------------
# Project bootstrap
# ---------------------------------------------------------------------------

# Initialise a fresh Asimov project.  Set PROJECT_NAME and ASIMOV_REGIME on
# the command line if needed (defaults: project directory name, production).
PROJECT_NAME  ?= $(notdir $(abspath $(PROJECT)))
ASIMOV_REGIME ?= production

$(PROJECT)/.asimov/ledger.yml:
	mkdir -p $(PROJECT)
	cd $(PROJECT); \
	$(ASIMOV) init "$(PROJECT_NAME)"; \
	$(ASIMOV) configuration update pipelines/environment $(ENV); \
	$(ASIMOV) configuration update general/webroot $(WEB_PATH); \
	$(ASIMOV) configuration update condor/user $(ACCOUNTING_USER); \
	$(ASIMOV) apply -f https://git.ligo.org/asimov/data/-/raw/main/defaults/production-pe.yaml; \
	$(ASIMOV) apply -f https://git.ligo.org/asimov/data/-/raw/main/defaults/cbcflow/$(ASIMOV_REGIME).yaml; \
	$(ASIMOV) apply -f https://git.ligo.org/asimov/data/-/raw/main/defaults/production-pe-priors.yaml

# Clone the CBCFlow metadata library.
# Requires CBCFLOW_TOKEN (personal access token) and CBCFLOW_GIT_URL to be set
# in the environment — never hardcode credentials in this file.
$(PROJECT)/cbcflow-library: $(PROJECT)/.asimov/ledger.yml
	@test -n "$(CBCFLOW_TOKEN)"  || (echo "ERROR: CBCFLOW_TOKEN is not set"; exit 1)
	@test -n "$(CBCFLOW_GIT_URL)" || (echo "ERROR: CBCFLOW_GIT_URL is not set"; exit 1)
	cd $(PROJECT); \
	git clone https://$(CBCFLOW_TOKEN)@$(CBCFLOW_GIT_URL) cbcflow-library; \
	cd cbcflow-library && ./setup-cbcflow-merge-strategy.sh

# ---------------------------------------------------------------------------
# Event management
# ---------------------------------------------------------------------------

# Add all events listed in event_list.txt:  make add-events
add-events: event_list.txt $(PROJECT)/.asimov/ledger.yml
	for event in $$(cat $<); do \
		$(MAKE) add-event-$$event; \
	done

# Add a single event and apply all blueprints from BLUEPRINTS_DIR:
#   make add-event-S230601bf
add-event-%: $(PROJECT)/.asimov/ledger.yml $(PROJECT)/cbcflow-library
	cd $(PROJECT); \
	$(ASIMOV) apply -p cbcflow -e $*; \
	for yaml in ../$(BLUEPRINTS_DIR)/*.yaml; do \
		$(ASIMOV) apply --event $* --file $$yaml; \
	done

# ---------------------------------------------------------------------------
# Stage 0: review combined summary pages (non-rejected, cheap settings)
# ---------------------------------------------------------------------------
combined_summary.mk: $(PROJECT)/.asimov/ledger.yml
	cd $(PROJECT)
	${PYTHON} ../scripts/generate_combining_command.py > ../$@

condor/review_summarypages.sub: $(PROJECT)/.asimov/ledger.yml
	mkdir -p condor
	cd $(PROJECT)
	${PYTHON} ../scripts/make_combining_submit_file.py > ../$@

# ---------------------------------------------------------------------------
# Stage 1: MVP SpinTaylor-only release (copy + relabel, no Condor needed)
# ---------------------------------------------------------------------------
spintaylor_release.mk: $(PROJECT)/.asimov/ledger.yml
	cd $(PROJECT)
	${PYTHON} ../scripts/generate_spintaylor_release_command.py > ../$@

# ---------------------------------------------------------------------------
# Stage 2A: release-quality summarypages (APPROVED only, spin evolution on)
# ---------------------------------------------------------------------------
release_summary.mk: $(PROJECT)/.asimov/ledger.yml
	cd $(PROJECT)
	${PYTHON} ../scripts/generate_release_summarypages.py > ../$@

condor/release_summarypages.sub: $(PROJECT)/.asimov/ledger.yml
	mkdir -p condor
	cd $(PROJECT)
	${PYTHON} ../scripts/make_release_summarypages_submit.py > ../$@

# ---------------------------------------------------------------------------
# Stage 2B: mix waveforms into combined release HDF5
# ---------------------------------------------------------------------------
mixing_release.mk: $(PROJECT)/.asimov/ledger.yml
	cd $(PROJECT)
	${PYTHON} ../scripts/generate_mixing_command.py > ../$@

condor/mixing_release.sub: $(PROJECT)/.asimov/ledger.yml
	mkdir -p condor
	cd $(PROJECT)
	${PYTHON} ../scripts/make_mixing_submit_file.py > ../$@

# ---------------------------------------------------------------------------
# Stage 2C: skymaps for full combined release (written directly to skymaps.mk)
# ---------------------------------------------------------------------------
skymaps.mk: $(PROJECT)/.asimov/ledger.yml
	cd $(PROJECT)
	${PYTHON} ../scripts/skymaps.py .asimov/ledger.yml > ../$@

# ---------------------------------------------------------------------------
# Convenience targets
# ---------------------------------------------------------------------------
makefiles: combined_summary.mk release_summary.mk spintaylor_release.mk \
           mixing_release.mk skymaps.mk

condor-files: condor/review_summarypages.sub \
              condor/release_summarypages.sub \
              condor/mixing_release.sub

# Quickview comparison plot for a single event: make quickview EVENT=S230601bf
quickview: ${WEB_PATH}/${EVENT}/quickview.png
${WEB_PATH}/%/quickview.png: $(PROJECT)/.asimov/ledger.yml
	cd $(PROJECT)
	${PYTHON} ../scripts/gwtc4-comp.py $*

tests:
	cd $(PROJECT)
	${PYTHON} -m unittest discover tests

# Verify every expected release file exists and is readable
release-check: $(PROJECT)/.asimov/ledger.yml
	cd $(PROJECT)
	${PYTHON} ../tests/test_release_files.py \
		--ledger .asimov/ledger.yml \
		--release-dir ../releases/$(RELEASE)

# Generate HTML + Markdown report summarising the release contents
release-report: $(PROJECT)/.asimov/ledger.yml
	cd $(PROJECT)
	${PYTHON} ../scripts/generate_release_report.py \
		--ledger .asimov/ledger.yml \
		--release-dir ../releases/$(RELEASE)

.PHONY: tests makefiles condor-files quickview release-check release-report add-events
