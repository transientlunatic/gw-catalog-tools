"""
Tests of the checks of the automated review.

These use synthetic data with known answers, and need neither asimov nor a
project.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from catalog_tools import review as utils


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
        assert result["details"]["single_detector_dominated"]
        assert not result["details"]["consistent_with_noise"]
        assert not result["details"]["implausibly_loud"]

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


class TestRailing:
    rng = np.random.default_rng(6)

    def test_a_posterior_which_fits_inside_its_prior_does_not_rail(self):
        low, high, _, bounds = utils.railing_check(self.rng.normal(50, 5, 20000))
        assert (low, high, bounds) == (False, False, (None, None))

    def test_a_posterior_cut_off_by_the_upper_bound_rails(self):
        samples = self.rng.normal(50, 20, 40000)
        samples = samples[samples < 60]
        low, high, _, (bound_low, bound_high) = utils.railing_check(samples)
        assert high and not low
        assert bound_low is None and bound_high > samples.max()

    def test_a_posterior_cut_off_by_the_lower_bound_rails(self):
        samples = self.rng.normal(50, 20, 40000)
        samples = samples[samples > 40]
        low, high, _, (bound_low, _) = utils.railing_check(samples)
        assert low and not high
        assert bound_low < samples.min()

    def test_mass_ratio_may_reach_one(self):
        # q = 1 is a genuine boundary, so piling up there is not railing
        samples = 1 - np.abs(self.rng.normal(0, 0.1, 20000))
        _, high, _, _ = utils.railing_check(samples, label="mass_ratio")
        assert not high
        _, high, _, _ = utils.railing_check(samples, label="chirp_mass")
        assert high

    def test_a_prior_which_already_extends_further_is_not_railing(self):
        samples = np.abs(self.rng.normal(0, 100, 40000)) + 10
        assert utils.railing_check(samples)[0]
        assert not utils.railing_check(samples, input_min_edge=1.0)[0]

    def test_agrees_with_peconfigurator(self):
        configurator = pytest.importorskip("peconfigurator.auxiliary.run_analyzer")
        for scale in (5, 20, 60):
            samples = self.rng.normal(50, scale, 40000)
            samples = samples[(samples > 20) & (samples < 80)]
            for kwargs in ({}, {"input_min_edge": 5.0}):
                ours = utils.railing_check(samples, **kwargs)
                theirs = configurator.railing_check(samples, Nbin=50, tolerance=2.0,
                                                    label="chirp_mass", **kwargs)
                assert ours[:2] == theirs[:2]
                assert ours[3][0] == theirs[3][0]

    def test_check_railing_reports_the_suggested_bounds(self):
        samples = self.rng.normal(50, 20, 40000)
        samples = samples[samples < 60]
        result = utils.check_railing(posterior(chirp_mass=samples), {}, "chirp_mass")
        assert result["status"] == "FAIL"
        assert result["details"]["suggested_upper_bound"] > 60
        ok = utils.check_railing(posterior(chirp_mass=self.rng.normal(50, 5, 20000)), {}, "chirp_mass")
        assert ok["status"] == "PASS"
        assert utils.check_railing(posterior(other=np.ones(10)), {}, "chirp_mass")["status"] == "SKIP"


def write_samples(path, **columns):
    """Write a minimal PESummary-style samples file."""
    import h5py

    n = len(next(iter(columns.values())))
    samples = np.zeros(n, dtype=[(name, float) for name in columns])
    for name, values in columns.items():
        samples[name] = values
    with h5py.File(path, "w") as handle:
        group = handle.create_group("analysis")
        group.create_dataset("posterior_samples", data=samples)
        group.create_group("priors/samples")
        group.create_dataset("meta_data/sampler/ln_bayes_factor", data=[12.0])
        group.create_dataset("meta_data/other/run_statistics/neffsamples", data=[5000.0])


class TestCommandLine:
    @pytest.fixture
    def good(self, tmp_path):
        rng = np.random.default_rng(7)
        n = 5000
        path = tmp_path / "good.h5"
        write_samples(path, chirp_mass=rng.normal(30, 3, n), luminosity_distance=rng.normal(900, 100, n),
                      network_matched_filter_snr=np.full(n, 8.0))
        return str(path)

    def test_review_passes_a_sensible_posterior(self, good, capsys):
        from catalog_tools import cli
        assert cli.main(["review", good]) == 0
        assert "PASS" in capsys.readouterr().out

    def test_review_json(self, good, capsys):
        import json
        from catalog_tools import cli
        cli.main(["review", good, "--json"])
        results = json.loads(capsys.readouterr().out)
        assert results["chirp_mass_railing"]["status"] == "PASS"
        assert results["effective_sample_size"]["details"]["neff"] == 5000.0

    def test_review_fails_a_railing_posterior(self, tmp_path, capsys):
        from catalog_tools import cli
        rng = np.random.default_rng(8)
        samples = rng.normal(50, 20, 40000)
        path = tmp_path / "railing.h5"
        write_samples(path, chirp_mass=samples[samples < 60])
        assert cli.main(["review", str(path)]) == 1
        assert "railing" in capsys.readouterr().out

    def test_review_reports_an_unreadable_file(self, tmp_path, capsys):
        from catalog_tools import cli
        assert cli.main(["review", str(tmp_path / "missing.h5")]) == 2

    def test_install_tests_writes_a_module_and_will_not_overwrite_it(self, tmp_path):
        from catalog_tools import cli
        directory = tmp_path / "tests"
        assert cli.main(["install-tests", str(directory)]) == 0
        assert "catalog_tools.asimov_tests" in (directory / "test_catalog_tools.py").read_text()
        assert cli.main(["install-tests", str(directory)]) == 1
        assert cli.main(["install-tests", str(directory), "--force"]) == 0


class TestFindingSamples:
    def test_candidates_are_in_order_of_preference(self):
        paths = utils.metafile_candidates("GW1", "run-2", webroot="/web", project_root="/proj",
                                          storage="results")
        assert paths == [
            "/web/GW1/run-2/pesummary/samples/posterior_samples.h5",
            "results/GW1/run-2/posterior_samples.h5",
            "/proj/results/GW1/run-2/posterior_samples.h5",
        ]

    def test_results_in_the_current_project_are_found_when_the_config_names_another_copy(
            self, tmp_path, monkeypatch):
        # The config says the project is at /elsewhere, but the analysis being
        # tested is in the project in the current directory.
        monkeypatch.chdir(tmp_path)
        samples = tmp_path / "results" / "GW1" / "run-2" / "posterior_samples.h5"
        samples.parent.mkdir(parents=True)
        samples.write_text("")
        monkeypatch.setattr(utils, "candidate_metafiles", lambda event, production: utils.metafile_candidates(
            event, production, webroot="/nonexistent/web", project_root="/elsewhere", storage="results"))
        assert utils.get_metafile_path("GW1", "run-2") == "results/GW1/run-2/posterior_samples.h5"

    def test_the_first_candidate_is_reported_when_nothing_exists(self, monkeypatch):
        monkeypatch.setattr(utils, "candidate_metafiles", lambda event, production: utils.metafile_candidates(
            event, production, webroot="/nonexistent/web", project_root="/elsewhere", storage="results"))
        assert utils.get_metafile_path("GW1", "run-2").startswith("/nonexistent/web/")
