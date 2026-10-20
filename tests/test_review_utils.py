"""
Tests of the helpers used by the automated review.

These use synthetic data with known answers, and need neither asimov nor a
project.
"""

from types import SimpleNamespace

import numpy as np
import pytest

import review_utils as utils


def production(pipeline="bilby", status="uploaded", review=None):
    """A stand-in for an asimov production."""
    return SimpleNamespace(
        pipeline=SimpleNamespace(name=pipeline),
        status_str=status,
        review=SimpleNamespace(status=review),
    )


def posterior(**columns):
    n = len(next(iter(columns.values())))
    array = np.zeros(n, dtype=[(name, float) for name in columns])
    for name, values in columns.items():
        array[name] = values
    return array


class TestSelection:
    @pytest.mark.parametrize("name", ["bilby", "Bilby", "BILBY"])
    def test_bilby_matches_any_case(self, name):
        # asimov 0.6 calls the pipeline "Bilby" and 0.7 calls it "bilby"
        assert utils.is_bilby(production(pipeline=name))

    def test_other_pipelines_do_not_match(self):
        assert not utils.is_bilby(production(pipeline="bayeswave"))

    @pytest.mark.parametrize("review, expected", [
        (None, True), ("APPROVED", False), ("REJECTED", False), ("DEPRECATED", False)])
    def test_only_unreviewed_analyses_are_reviewed(self, review, expected, monkeypatch):
        monkeypatch.delenv("PE_REVIEW_INCLUDE_REVIEWED", raising=False)
        assert utils.should_review(production(review=review)) is expected

    def test_reviewed_analyses_can_be_included(self, monkeypatch):
        monkeypatch.setenv("PE_REVIEW_INCLUDE_REVIEWED", "1")
        assert utils.should_review(production(review="REJECTED"))

    def test_analyses_which_have_not_been_uploaded_are_not_reviewed(self):
        assert not utils.should_review(production(status="running"))

    def test_approval_is_case_insensitive(self):
        assert utils.is_approved(production(review="approved"))
        assert not utils.is_approved(production(review=None))


class TestModes:
    rng = np.random.default_rng(0)

    def test_gaussians_have_one_mode(self):
        for _ in range(10):
            assert len(utils.find_modes(self.rng.normal(size=15000))) == 1

    def test_skewed_and_bounded_posteriors_have_one_mode(self):
        for samples in (self.rng.lognormal(0, 0.8, 15000), self.rng.gamma(2, size=15000),
                        np.abs(self.rng.normal(size=15000)), self.rng.beta(2, 0.7, 15000)):
            assert len(utils.find_modes(samples)) == 1

    def test_well_separated_modes_are_found(self):
        samples = np.r_[self.rng.normal(-6, 1, 7500), self.rng.normal(6, 1, 7500)]
        modes = utils.find_modes(samples)
        assert len(modes) == 2
        assert modes[0]["location"] == pytest.approx(-6, abs=0.5)
        assert modes[1]["location"] == pytest.approx(6, abs=0.5)
        assert modes[0]["mass"] == pytest.approx(0.5, abs=0.05)

    def test_modes_with_no_gap_between_them_are_found(self):
        # The modes at 25 and 30 overlap, so there is no gap in the samples,
        # which the earlier gap-based test needed.
        samples = np.r_[self.rng.normal(25, 1.5, 7500), self.rng.normal(30, 1.5, 7500)]
        assert len(utils.find_modes(samples)) == 2

    def test_a_minor_mode_is_ignored(self):
        samples = np.r_[self.rng.normal(0, 1, 14700), self.rng.normal(5, 1, 300)]
        assert len(utils.find_modes(samples)) == 1

    def test_log_space_for_scale_like_parameters(self):
        samples = np.exp(np.r_[self.rng.normal(np.log(100), 0.1, 7500),
                               self.rng.normal(np.log(1000), 0.1, 7500)])
        modes = utils.find_modes(samples, log=True)
        assert [m["location"] for m in modes] == pytest.approx([100, 1000], rel=0.05)

    def test_too_few_samples_or_no_spread(self):
        assert utils.find_modes(np.ones(1000)) == []
        assert utils.find_modes(np.arange(10.0)) == []


class TestSkyArea:
    def test_whole_sky_is_about_41000_square_degrees(self):
        rng = np.random.default_rng(1)
        ra = rng.uniform(0, 2 * np.pi, 200000)
        dec = np.arcsin(rng.uniform(-1, 1, 200000))
        area, cell = utils.calculate_sky_area_90(ra, dec)
        assert area == pytest.approx(0.9 * 41253, rel=0.03)
        assert cell == pytest.approx(41253 / (64 * 32), rel=0.01)

    def test_a_compact_patch_is_small(self):
        rng = np.random.default_rng(2)
        ra = rng.normal(1.0, 0.05, 20000)
        dec = rng.normal(0.3, 0.05, 20000)
        area, _ = utils.calculate_sky_area_90(ra, dec)
        assert area < 300


class TestChecks:
    prior = {"a_1": np.random.default_rng(3).uniform(0, 0.99, 20000)}

    def test_spins_pinned_to_the_bound_are_flagged(self):
        post = posterior(a_1=np.full(1000, 0.98))
        assert utils.check_spin_prior(post, self.prior)["status"] == "WARNING"

    def test_spins_which_follow_the_prior_pass(self):
        post = posterior(a_1=np.random.default_rng(4).uniform(0, 0.99, 5000))
        assert utils.check_spin_prior(post, self.prior)["status"] == "PASS"

    def test_effective_sample_size_prefers_the_sampler_estimate(self):
        result = utils.check_effective_sample_size(20000, {"neffsamples": 500})
        assert result["status"] == "FAIL"
        assert utils.check_effective_sample_size(500, {})["status"] == "FAIL"
        assert utils.check_effective_sample_size(20000, {"neffsamples": 8000})["status"] == "PASS"

    @pytest.mark.parametrize("snr, ln_bf, flagged", [
        (7.5, 10.0, False),     # a typical near-threshold candidate
        (0.8, -0.08, True),     # consistent with noise
        (32.0, 450.0, True),    # far too loud
    ])
    def test_snr_plausibility(self, snr, ln_bf, flagged):
        post = posterior(network_matched_filter_snr=np.full(100, snr))
        result = utils.check_snr_plausibility(post, {"ln_bayes_factor": ln_bf})
        assert (result["status"] == "WARNING") is flagged

    def test_a_single_dominant_detector_is_flagged(self):
        post = posterior(network_matched_filter_snr=np.full(100, 9.0),
                         H1_optimal_snr=np.full(100, 3.0), L1_optimal_snr=np.full(100, 9.0))
        result = utils.check_snr_plausibility(post, {"ln_bayes_factor": 10.0})
        assert result["status"] == "WARNING" and "L1" in result["message"]

    def test_inconsistent_derived_parameters_fail(self):
        m1, m2 = np.full(100, 30.0), np.full(100, 20.0)
        good = posterior(mass_1=m1, mass_2=m2, total_mass=m1 + m2)
        assert utils.check_derived_parameters(good)["status"] == "PASS"
        bad = posterior(mass_1=m1, mass_2=m2, total_mass=m1 + m2 + 1)
        assert utils.check_derived_parameters(bad)["status"] == "FAIL"

    def test_network_snr_must_be_the_quadrature_sum(self):
        good = posterior(network_optimal_snr=np.full(10, 5.0),
                         H1_optimal_snr=np.full(10, 3.0), L1_optimal_snr=np.full(10, 4.0))
        assert utils.check_network_snr(good)["status"] == "PASS"
        bad = posterior(network_optimal_snr=np.full(10, 6.0),
                        H1_optimal_snr=np.full(10, 3.0), L1_optimal_snr=np.full(10, 4.0))
        assert utils.check_network_snr(bad)["status"] == "FAIL"

    def test_a_nearly_constant_likelihood_fails(self):
        assert utils.check_likelihood(posterior(log_likelihood=np.full(100, 5.0)))["status"] == "FAIL"
        varied = np.random.default_rng(5).normal(20, 2, 100)
        assert utils.check_likelihood(posterior(log_likelihood=varied))["status"] == "PASS"
