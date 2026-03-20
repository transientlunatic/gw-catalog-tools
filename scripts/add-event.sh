event=$1
#asimov apply -f ../blueprints/events/gwtc-4/$event.yaml
# asimov apply \
#        --event $event \
#        --file ../blueprints/analyses/psd-from-release.yaml
# asimov apply \
#        --event $event \
#        --file ../blueprints/analyses/bilby-bbh/analysis_bilby_IMRPhenomXPHM-SpinTaylor.yaml
# asimov apply \
#        --event $event \
#        --file ../blueprints/analyses/bilby-bbh/analysis_bilby_SEOBNRv5PHM.yaml
# asimov apply \
#        --event $event \
#        --file ../blueprints/analyses/bilby-bbh/analysis_bilby_NRSur7Dq4.yaml
asimov apply \
       --event $event \
       --file ../blueprints/analyses/bilby-bbh/analysis_bilby_IMRPhenomXO4a.yaml
