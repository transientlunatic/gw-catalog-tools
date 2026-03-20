from asimov.testing import AsimovTest
import os
import sys
import peconfigurator.auxiliary.run_analyzer as configurator

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from ledger import open_pesummary_samples


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
