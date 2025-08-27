from asimov.testing import AsimovTest
import configparser
from bilby_pipe import bilbyargparser
import ast
import os
import glob
import dill
import h5py as h5
import numpy as np
import peconfigurator.auxiliary.run_analyzer as configurator

def open_pesummary_samples(metafile):
    with h5.File(metafile) as metafile_h:
        analyses = list(metafile_h.keys())
        if "history" in analyses:
            analyses.remove("history")
        if "version" in analyses:
            analyses.remove("version")

        key = analyses[0]

        posterior = np.array(metafile_h[key]['posterior_samples'])
        prior = dict(metafile_h[key]['priors']['samples'])
        prior = {name: np.array(data) for name, data in prior.items()}
        #print(prior)
        
    return posterior, prior


class TestBilbyAnalysesForRailing(AsimovTest):
    def test_distance_railing(self):
        """Check for railing in the distance posterior"""

        for event in self.events:
            for production in event.productions:
                if (production.pipeline.name.lower() == "bilby") \
                   and (production.status_str.lower() == "uploaded")\
                   and (production.review.status not in ["REJECTED", "APPROVED", "DEPRECATED"]):
                    with self.subTest(event=event.name, production=production.name):
                        
                        metafile = f"/home/pe.o4/public_html/gwtc4-fogg/{event.name}/{production.name}/pesummary/samples/posterior_samples.h5"
                        posterior, prior = open_pesummary_samples(metafile)
                        railing_down, railing_up, rec, _ = configurator.railing_check(posterior["luminosity_distance"],
                                                             Nbin=50,
                                                             tolerance=2.0,
                                                             label="luminosity_distance")
                        self.assertFalse(railing_up, msg=f"Railing detected in luminosity distance at the upper bound; recommend this is adjusted to {_[1]}")
                        #self.assertFalse(railing_down)

    def test_chirp_railing(self):
        """Check for railing in the chirp_mass posterior"""

        for event in self.events:
            for production in event.productions:
                if (production.pipeline.name.lower() == "bilby") \
                   and (production.status_str.lower() == "uploaded") \
                   and (production.review.status not in ["REJECTED", "APPROVED", "DEPRECATED"]):
                    with self.subTest(event=event.name, production=production.name):
                        
                        metafile = f"/home/pe.o4/public_html/gwtc4-fogg/{event.name}/{production.name}/pesummary/samples/posterior_samples.h5"
                        posterior, prior = open_pesummary_samples(metafile)
                        railing_down, railing_up, rec, _ = configurator.railing_check(posterior["chirp_mass"],
                                                             Nbin=50,
                                                             tolerance=2.0,
                                                             label="chirp_mass")
                        self.assertFalse(railing_up, msg=f"Railing detected in mass ratio at the upper bound; recommend this is adjusted to {_[1]}")
                        self.assertFalse(railing_down, msg=f"Railing detected in mass ratio at the lower bound; recommend this is adjusted to {_[0]}")
