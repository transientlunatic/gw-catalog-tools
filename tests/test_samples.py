from asimov.testing import AsimovTest
from bilby_pipe import bilbyargparser
import glob
import dill

class TestBilbyAnalyses(AsimovTest):
    def test_sample_loss(self):
        """Verify the component chirp mass prior is correct."""

        for event in self.events:
            for production in event.productions:
                if (production.pipeline.name == "Bilby") \
                   and (production.status_str.lower() == "uploaded") \
                   and (production.review.status == "APPROVED"):
                    with self.subTest(event=event.name, production=production.name):
                        # Read the config file
                        pickles = glob.glob(f"{production.rundir}/result/*dynesty.pickle")

                        for fname in pickles:
                            with open(fname, "rb") as ff:
                                result = dill.load(ff)
                                self.assertEqual((result.niter + result.nlive), len(result.samples))

    def test_sample_loss_unreviewed(self):
        """Verify samples are not lost in unreviewed analyses"""

        for event in self.events:
            for production in event.productions:
                if (production.pipeline.name == "Bilby") \
                   and (production.status_str.lower() == "uploaded") \
                   and len(production.review) == 0:
                    with self.subTest(event=event.name, production=production.name):
                        # Read the config file
                        pickles = glob.glob(f"{production.rundir}/result/*dynesty.pickle")

                        for fname in pickles:
                            with open(fname, "rb") as ff:
                                result = dill.load(ff)
                                self.assertEqual((result.niter + result.nlive), len(result.samples))

                                

    def test_sample_loss_reviewed(self):
        """Verify samples are not lost in reviewed analyses"""

        for event in self.events:
            for production in event.productions:
                if (production.pipeline.name == "Bilby") \
                   and (production.status_str.lower() == "uploaded") \
                   and production.review.status is not None and production.review.status.lower() == "approved":
                    with self.subTest(event=event.name, production=production.name):
                        # Read the config file
                        pickles = glob.glob(f"{production.rundir}/result/*dynesty.pickle")

                        for fname in pickles:
                            with open(fname, "rb") as ff:
                                result = dill.load(ff)
                                self.assertEqual((result.niter + result.nlive), len(result.samples))

                                
