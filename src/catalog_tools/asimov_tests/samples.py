"""Check that no samples were lost in the sampler output."""

import glob

import dill
from asimov.testing import AsimovTest

from catalog_tools.review import is_approved, is_bilby, is_uploaded, should_review


class TestSampleLoss(AsimovTest):
    def sample_loss(self, selector):
        """
        Verify that the dynesty result of each selected analysis holds all of
        the samples: ``niter + nlive`` of them.

        Failures are reported together with their event and analysis. The test
        is skipped if no result was found, so that a test which has checked
        nothing can't pass.
        """
        checked = 0
        failures = []
        for event in self.events:
            for production in event.productions:
                if not (is_bilby(production) and is_uploaded(production) and selector(production)):
                    continue
                for fname in glob.glob(f"{production.rundir}/result/*dynesty.pickle"):
                    with open(fname, "rb") as ff:
                        result = dill.load(ff)
                    checked += 1
                    expected = result.niter + result.nlive
                    if expected != len(result.samples):
                        failures.append(f"{event.name} {production.name}: {len(result.samples)} samples, "
                                        f"expected niter + nlive = {expected}")
        if failures:
            self.fail("\n".join(failures))
        if not checked:
            self.skipTest("No sampler results found for the selected analyses")

    def test_sample_loss_approved(self):
        """Verify samples are not lost in approved analyses."""
        self.sample_loss(is_approved)

    def test_sample_loss_unreviewed(self):
        """Verify samples are not lost in analyses which are awaiting review."""
        self.sample_loss(should_review)
