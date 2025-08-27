
.ONESHELL:
SHELL=/usr/bin/bash
ENV=/home/pe.o4/GWTC4-fogg/environment
WEB_PATH= /home/pe.o4/public_html/gwtc4

SUMMARYPAGES=${ENV}/bin/summarypages
SUMMARY_COMBINE_POSTERIORS=${ENV}/bin/summarycombine_posteriors
SUMMARYMODIFY=${ENV}/bin/summarymodify
EVENTS=$(shell ls working)

WEB_PATH= /home/pe.o4/public_html/gwtc4-fogg
RELEASE=test
ENV=/home/pe.o4/GWTC4-fogg/environment
SUMMARYPAGES=${ENV}/bin/summarypages
PYTHON=${ENV}/bin/python

RELEASE=staging

include combined_summary.mk

# Make the combined summarypages makefile
combined_summary.mk: project/.asimov/ledger.yml
	cd project
	${PYTHON} ../scripts/generate_combining_command.py > ../$@

# Make a submit file for all the combined pages
condor/all_combined.sub: combined_summary.mk
	mkdir -p condor
	cd project
	${PYTHON} ../scripts/make_combining_submit_file.py > ../$@

${WEB_PATH}/%/quickview.png: project/.asimov/ledger.yml
	cd project
	${PYTHON} ../scripts/gwtc4-comp.py $*


tests:
	cd project
	${PYTHON} -m unittest discover tests

makefiles: combined_summary.mk

.PHONY: tests makefiles
