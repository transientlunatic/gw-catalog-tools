"""Check for railing of the posteriors against the prior bounds."""

import numpy as np
from asimov.testing import AsimovTest
import peconfigurator.auxiliary.run_analyzer as configurator

import os

from review_utils import get_metafile_path, open_pesummary_samples, should_review


class TestBilbyAnalysesForRailing(AsimovTest):
    """
    A posterior which piles up against a prior bound means the bound is cutting
    off real support, so the prior should be widened.  The check is the one used
    by peconfigurator: the end bins of a histogram are compared to the highest.
    """

    def railing_failures(self, param, check_lower, check_upper, prior_minimum=False):
        """
        Check ``param`` in every analysis which should be reviewed, and fail if
        the posterior rails at an end which is being checked.

        Every analysis is checked and the failures reported together, each
        with its event and analysis.  If ``prior_minimum`` is set the smallest
        prior sample is passed to the check, which then ignores lower railing
        when the prior already extends below the suggested new bound.
        """
        failures = []
        checked = 0
        for event in self.events:
            for production in event.productions:
                if not should_review(production):
                    continue
                label = f"{event.name} {production.name}"
                metafile = get_metafile_path(event.name, production.name)
                if not os.path.exists(metafile):
                    failures.append(f"{label}: samples not found at {metafile}")
                    continue

                posterior, prior = open_pesummary_samples(metafile)
                kwargs = {}
                if prior_minimum and len(prior.get(param, [])):
                    kwargs["input_min_edge"] = float(np.min(prior[param]))

                down, up, _, (bound_low, bound_high) = configurator.railing_check(
                    posterior[param], Nbin=50, tolerance=2.0, label=param, **kwargs)
                checked += 1

                if check_lower and down:
                    failures.append(f"{label}: railing in {param} at the lower bound; "
                                    f"recommend this is adjusted to {bound_low}")
                if check_upper and up:
                    failures.append(f"{label}: railing in {param} at the upper bound; "
                                    f"recommend this is adjusted to {bound_high}")
        if not checked and not failures:
            self.skipTest("No analyses to review")
        if failures:
            self.fail("\n".join(failures))

    def test_distance_railing(self):
        """Check for railing in the luminosity distance posterior."""
        # Distance can legitimately pile up towards zero, which the check
        # ignores if the prior already extends below the suggested bound.
        self.railing_failures("luminosity_distance", check_lower=True, check_upper=True,
                              prior_minimum=True)

    def test_chirp_railing(self):
        """Check for railing in the chirp mass posterior."""
        self.railing_failures("chirp_mass", check_lower=True, check_upper=True)
