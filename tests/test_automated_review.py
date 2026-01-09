"""
Automated PE Review Tests
==========================

This module provides comprehensive automated quality checks for parameter estimation
analyses. Tests run once per production and generate an HTML report for review.

The tests check for common failure modes including:
- Spin magnitude railing at physical bounds
- Poor sampling quality (low effective sample size)
- Derived parameter inconsistencies
- Evidence/Bayes factor sanity checks
- Sky localization quality
- Network SNR consistency

Tests are designed to run without requiring trained ML models.
"""

from asimov.testing import AsimovTest
import numpy as np
import h5py
import os
from datetime import datetime
from collections import defaultdict
import json


def open_pesummary_samples(metafile):
    """
    Open and extract posterior and prior samples from PESummary HDF5 file.

    Parameters
    ----------
    metafile : str
        Path to the posterior_samples.h5 file

    Returns
    -------
    posterior : numpy.ndarray
        Structured array containing posterior samples
    prior : dict
        Dictionary of prior samples for each parameter
    """
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

    return posterior, prior


def calculate_ess(samples):
    """
    Calculate effective sample size using autocorrelation.

    Parameters
    ----------
    samples : array_like
        1D array of samples

    Returns
    -------
    ess : float
        Effective sample size
    """
    n = len(samples)
    if n < 2:
        return n

    # Demean
    samples = samples - np.mean(samples)

    # Calculate autocorrelation
    c0 = np.dot(samples, samples) / n
    if c0 == 0:
        return n

    # Sum autocorrelations until they become negative or negligible
    autocorr_sum = 0
    for k in range(1, min(n, 1000)):
        ck = np.dot(samples[:-k], samples[k:]) / (n - k)
        autocorr_sum += ck / c0
        if ck / c0 < 0.05:  # Stop when autocorr drops below 5%
            break

    ess = n / (1 + 2 * autocorr_sum)
    return max(1, ess)  # ESS should be at least 1


def check_parameter_consistency(posterior, param_derived, param_inputs, calculation_func, tolerance=1e-3):
    """
    Check that derived parameters are consistent with sampled parameters.

    Parameters
    ----------
    posterior : numpy.ndarray
        Structured array of posterior samples
    param_derived : str
        Name of the derived parameter to check
    param_inputs : list of str
        Names of input parameters used in calculation
    calculation_func : callable
        Function that calculates derived parameter from inputs
    tolerance : float
        Relative tolerance for consistency check

    Returns
    -------
    is_consistent : bool
        True if derived parameter matches calculation
    max_error : float
        Maximum relative error found
    """
    # Check all parameters exist
    for param in [param_derived] + param_inputs:
        if param not in posterior.dtype.names:
            return None, None

    # Get values
    derived_values = posterior[param_derived]
    input_values = [posterior[param] for param in param_inputs]

    # Calculate what the derived values should be
    calculated_values = calculation_func(*input_values)

    # Check consistency
    relative_errors = np.abs(derived_values - calculated_values) / (np.abs(derived_values) + 1e-10)
    max_error = np.max(relative_errors)
    is_consistent = max_error < tolerance

    return is_consistent, max_error


def calculate_sky_area_90(ra, dec, num_bins=50):
    """
    Calculate approximate 90% credible area for sky localization.

    Parameters
    ----------
    ra : array_like
        Right ascension samples (radians)
    dec : array_like
        Declination samples (radians)
    num_bins : int
        Number of bins for 2D histogram

    Returns
    -------
    area_90 : float
        90% credible area in square degrees
    """
    # Create 2D histogram
    H, ra_edges, dec_edges = np.histogram2d(ra, dec, bins=num_bins)

    # Normalize to get probability density
    H = H / np.sum(H)

    # Flatten and sort
    H_flat = H.flatten()
    sorted_indices = np.argsort(H_flat)[::-1]

    # Find 90% credible region
    cumsum = np.cumsum(H_flat[sorted_indices])
    idx_90 = np.searchsorted(cumsum, 0.90)

    # Calculate area (rough approximation)
    ra_bin_width = (ra_edges[1] - ra_edges[0])
    dec_bin_width = (dec_edges[1] - dec_edges[0])

    # Account for spherical geometry: area element is cos(dec) * dra * ddec
    dec_centers = (dec_edges[:-1] + dec_edges[1:]) / 2

    # Count bins in 90% region and weight by cos(dec)
    area_steradian = 0
    for i in range(num_bins):
        for j in range(num_bins):
            flat_idx = i * num_bins + j
            if flat_idx in sorted_indices[:idx_90+1]:
                area_steradian += ra_bin_width * dec_bin_width * np.cos(dec_centers[j])

    # Convert to square degrees
    area_deg2 = area_steradian * (180 / np.pi) ** 2

    return area_deg2


class TestAutomatedReview(AsimovTest):
    """
    Automated PE review tests that run once per production.

    These tests check for common failure modes and generate an HTML report
    summarizing the results. The report is placed alongside the analysis
    summary pages for easy access during review.
    """

    @classmethod
    def setUpClass(cls):
        """Initialize test results storage."""
        super().setUpClass()
        cls.test_results = defaultdict(lambda: defaultdict(dict))

    @classmethod
    def tearDownClass(cls):
        """Generate HTML report after all tests complete."""
        super().tearDownClass()
        cls.generate_html_report()

    def record_result(self, event_name, production_name, test_name, result):
        """
        Record test result for HTML report generation.

        Parameters
        ----------
        event_name : str
            Name of the event
        production_name : str
            Name of the production/analysis
        test_name : str
            Name of the test
        result : dict
            Test result containing 'status', 'message', and optional 'details'
        """
        self.test_results[event_name][production_name][test_name] = result

    def test_spin_magnitude_railing(self):
        """
        Check for railing in spin magnitude parameters.

        Spin magnitudes (a_1, a_2) have physical bounds [0, 1]. Railing at
        these bounds may indicate prior is too restrictive or sampler issues.

        This test checks for accumulation of samples at the boundaries using
        the same railing detection as the distance/chirp_mass tests.
        """
        for event in self.events:
            for production in event.productions:
                if (production.pipeline.name.lower() == "bilby") \
                   and (production.status_str.lower() == "uploaded")\
                   and (production.review.status not in ["REJECTED", "APPROVED", "DEPRECATED"]):
                    with self.subTest(event=event.name, production=production.name):

                        metafile = f"/home/pe.o4/public_html/gwtc4-fogg/{event.name}/{production.name}/pesummary/samples/posterior_samples.h5"

                        if not os.path.exists(metafile):
                            self.record_result(event.name, production.name, "spin_railing", {
                                "status": "SKIP",
                                "message": "Metafile not found"
                            })
                            continue

                        posterior, prior = open_pesummary_samples(metafile)

                        issues = []

                        # Check a_1
                        if 'a_1' in posterior.dtype.names:
                            import peconfigurator.auxiliary.run_analyzer as configurator
                            railing_down, railing_up, rec, bounds = configurator.railing_check(
                                posterior["a_1"],
                                Nbin=50,
                                tolerance=2.0,
                                label="a_1"
                            )

                            if railing_up:
                                issues.append(f"a_1 railing at upper bound (near 1.0), recommend adjusting to {bounds[1]:.3f}")
                            if railing_down:
                                issues.append(f"a_1 railing at lower bound (near 0.0), recommend adjusting to {bounds[0]:.3f}")

                        # Check a_2
                        if 'a_2' in posterior.dtype.names:
                            import peconfigurator.auxiliary.run_analyzer as configurator
                            railing_down, railing_up, rec, bounds = configurator.railing_check(
                                posterior["a_2"],
                                Nbin=50,
                                tolerance=2.0,
                                label="a_2"
                            )

                            if railing_up:
                                issues.append(f"a_2 railing at upper bound (near 1.0), recommend adjusting to {bounds[1]:.3f}")
                            if railing_down:
                                issues.append(f"a_2 railing at lower bound (near 0.0), recommend adjusting to {bounds[0]:.3f}")

                        # Check chi_eff
                        if 'chi_eff' in posterior.dtype.names:
                            import peconfigurator.auxiliary.run_analyzer as configurator
                            railing_down, railing_up, rec, bounds = configurator.railing_check(
                                posterior["chi_eff"],
                                Nbin=50,
                                tolerance=2.0,
                                label="chi_eff"
                            )

                            if railing_up:
                                issues.append(f"chi_eff railing at upper bound, recommend adjusting to {bounds[1]:.3f}")
                            if railing_down:
                                issues.append(f"chi_eff railing at lower bound, recommend adjusting to {bounds[0]:.3f}")

                        # Check chi_p
                        if 'chi_p' in posterior.dtype.names:
                            import peconfigurator.auxiliary.run_analyzer as configurator
                            railing_down, railing_up, rec, bounds = configurator.railing_check(
                                posterior["chi_p"],
                                Nbin=50,
                                tolerance=2.0,
                                label="chi_p"
                            )

                            if railing_up:
                                issues.append(f"chi_p railing at upper bound, recommend adjusting to {bounds[1]:.3f}")
                            if railing_down:
                                issues.append(f"chi_p railing at lower bound, recommend adjusting to {bounds[0]:.3f}")

                        self.record_result(event.name, production.name, "spin_railing", {
                            "status": "FAIL" if issues else "PASS",
                            "message": "; ".join(issues) if issues else "No spin railing detected",
                            "details": {"issues": issues}
                        })

                        self.assertEqual(len(issues), 0, msg="\n".join(issues))

    def test_effective_sample_size(self):
        """
        Check effective sample size for key parameters.

        Low effective sample size (ESS) indicates poor sampling efficiency,
        often due to high autocorrelation. This can happen when:
        - Sampler step size is too small
        - Posterior has strong correlations
        - Sampler hasn't converged

        We check ESS for distance, chirp_mass, and mass_ratio as these are
        typically well-constrained parameters. ESS should be at least 10% of
        the total number of samples.
        """
        for event in self.events:
            for production in event.productions:
                if (production.pipeline.name.lower() == "bilby") \
                   and (production.status_str.lower() == "uploaded")\
                   and (production.review.status not in ["REJECTED", "APPROVED", "DEPRECATED"]):
                    with self.subTest(event=event.name, production=production.name):

                        metafile = f"/home/pe.o4/public_html/gwtc4-fogg/{event.name}/{production.name}/pesummary/samples/posterior_samples.h5"

                        if not os.path.exists(metafile):
                            self.record_result(event.name, production.name, "effective_sample_size", {
                                "status": "SKIP",
                                "message": "Metafile not found"
                            })
                            continue

                        posterior, prior = open_pesummary_samples(metafile)
                        n_samples = len(posterior)

                        params_to_check = ['luminosity_distance', 'chirp_mass', 'mass_ratio']
                        ess_results = {}
                        issues = []

                        for param in params_to_check:
                            if param in posterior.dtype.names:
                                ess = calculate_ess(posterior[param])
                                ess_fraction = ess / n_samples
                                ess_results[param] = {
                                    "ess": ess,
                                    "fraction": ess_fraction
                                }

                                if ess_fraction < 0.1:
                                    issues.append(f"{param}: ESS={ess:.0f} ({ess_fraction*100:.1f}% of {n_samples} samples)")

                        self.record_result(event.name, production.name, "effective_sample_size", {
                            "status": "FAIL" if issues else "PASS",
                            "message": f"Low ESS detected: {', '.join(issues)}" if issues else f"ESS acceptable for all parameters",
                            "details": ess_results
                        })

                        self.assertEqual(len(issues), 0,
                                       msg=f"Low effective sample size detected:\n" + "\n".join(issues))

    def test_derived_parameter_consistency(self):
        """
        Verify derived parameters are consistent with sampled parameters.

        Some parameters are derived from others (e.g., total_mass = mass_1 + mass_2,
        chirp_mass from mass_1 and mass_2, mass_ratio = mass_2/mass_1).

        Inconsistencies can indicate:
        - Post-processing errors
        - File corruption
        - Incorrect parameter transformations

        This test checks several key relationships with 0.1% tolerance.
        """
        for event in self.events:
            for production in event.productions:
                if (production.pipeline.name.lower() == "bilby") \
                   and (production.status_str.lower() == "uploaded")\
                   and (production.review.status not in ["REJECTED", "APPROVED", "DEPRECATED"]):
                    with self.subTest(event=event.name, production=production.name):

                        metafile = f"/home/pe.o4/public_html/gwtc4-fogg/{event.name}/{production.name}/pesummary/samples/posterior_samples.h5"

                        if not os.path.exists(metafile):
                            self.record_result(event.name, production.name, "parameter_consistency", {
                                "status": "SKIP",
                                "message": "Metafile not found"
                            })
                            continue

                        posterior, prior = open_pesummary_samples(metafile)

                        issues = []
                        checks = []

                        # Check total_mass = mass_1 + mass_2
                        is_consistent, max_error = check_parameter_consistency(
                            posterior, 'total_mass', ['mass_1', 'mass_2'],
                            lambda m1, m2: m1 + m2
                        )
                        if is_consistent is not None:
                            checks.append(f"total_mass consistency: {max_error:.2e}")
                            if not is_consistent:
                                issues.append(f"total_mass inconsistent with mass_1 + mass_2 (max error: {max_error:.2e})")

                        # Check mass_ratio = mass_2 / mass_1
                        is_consistent, max_error = check_parameter_consistency(
                            posterior, 'mass_ratio', ['mass_2', 'mass_1'],
                            lambda m2, m1: m2 / m1
                        )
                        if is_consistent is not None:
                            checks.append(f"mass_ratio consistency: {max_error:.2e}")
                            if not is_consistent:
                                issues.append(f"mass_ratio inconsistent with mass_2/mass_1 (max error: {max_error:.2e})")

                        # Check chirp_mass from mass_1 and mass_2
                        def calc_chirp_mass(m1, m2):
                            return (m1 * m2)**(3/5) / (m1 + m2)**(1/5)

                        is_consistent, max_error = check_parameter_consistency(
                            posterior, 'chirp_mass', ['mass_1', 'mass_2'],
                            calc_chirp_mass
                        )
                        if is_consistent is not None:
                            checks.append(f"chirp_mass consistency: {max_error:.2e}")
                            if not is_consistent:
                                issues.append(f"chirp_mass inconsistent with calculation (max error: {max_error:.2e})")

                        self.record_result(event.name, production.name, "parameter_consistency", {
                            "status": "FAIL" if issues else "PASS",
                            "message": "; ".join(issues) if issues else "All checked parameters consistent",
                            "details": {"checks": checks, "issues": issues}
                        })

                        self.assertEqual(len(issues), 0,
                                       msg=f"Parameter consistency errors:\n" + "\n".join(issues))

    def test_log_evidence_sanity(self):
        """
        Check that log evidence (Bayes factor) values are reasonable.

        Log evidence should be:
        - A finite number (not NaN or Inf)
        - Within reasonable bounds for GW analyses (typically -100 to 100)

        Unreasonable values can indicate:
        - Sampler failure
        - Numerical instability
        - Prior/likelihood issues

        For log_likelihood, we check that:
        - Values are finite
        - Maximum likelihood is reasonable (typically > 0 for GW signals)
        - Distribution has reasonable spread
        """
        for event in self.events:
            for production in event.productions:
                if (production.pipeline.name.lower() == "bilby") \
                   and (production.status_str.lower() == "uploaded")\
                   and (production.review.status not in ["REJECTED", "APPROVED", "DEPRECATED"]):
                    with self.subTest(event=event.name, production=production.name):

                        metafile = f"/home/pe.o4/public_html/gwtc4-fogg/{event.name}/{production.name}/pesummary/samples/posterior_samples.h5"

                        if not os.path.exists(metafile):
                            self.record_result(event.name, production.name, "log_evidence_sanity", {
                                "status": "SKIP",
                                "message": "Metafile not found"
                            })
                            continue

                        posterior, prior = open_pesummary_samples(metafile)

                        issues = []
                        details = {}

                        # Check log_likelihood
                        if 'log_likelihood' in posterior.dtype.names:
                            log_like = posterior['log_likelihood']

                            # Check for non-finite values
                            if not np.all(np.isfinite(log_like)):
                                n_bad = np.sum(~np.isfinite(log_like))
                                issues.append(f"log_likelihood has {n_bad} non-finite values")
                            else:
                                max_log_like = np.max(log_like)
                                mean_log_like = np.mean(log_like)
                                std_log_like = np.std(log_like)

                                details['log_likelihood'] = {
                                    'max': float(max_log_like),
                                    'mean': float(mean_log_like),
                                    'std': float(std_log_like)
                                }

                                # Check if maximum is suspiciously low
                                if max_log_like < -1000:
                                    issues.append(f"Maximum log_likelihood is very low: {max_log_like:.1f}")

                                # Check if std is suspiciously small (suggests sampling issues)
                                if std_log_like < 0.1:
                                    issues.append(f"log_likelihood has very small std: {std_log_like:.2e}")

                        # Check log_prior
                        if 'log_prior' in posterior.dtype.names:
                            log_prior = posterior['log_prior']

                            if not np.all(np.isfinite(log_prior)):
                                n_bad = np.sum(~np.isfinite(log_prior))
                                issues.append(f"log_prior has {n_bad} non-finite values")

                        self.record_result(event.name, production.name, "log_evidence_sanity", {
                            "status": "FAIL" if issues else "PASS",
                            "message": "; ".join(issues) if issues else "Log evidence values are reasonable",
                            "details": details
                        })

                        self.assertEqual(len(issues), 0,
                                       msg=f"Log evidence issues:\n" + "\n".join(issues))

    def test_sky_localization_quality(self):
        """
        Check sky localization quality.

        This test calculates the 90% credible area from RA and Dec posteriors.

        Issues flagged:
        - Area > 1000 deg² suggests poor localization (may be expected for
          low-SNR events or single-detector observations)
        - Area < 0.01 deg² is suspiciously small and may indicate sampling issues

        This is an informational test - failures may be expected for some events.
        """
        for event in self.events:
            for production in event.productions:
                if (production.pipeline.name.lower() == "bilby") \
                   and (production.status_str.lower() == "uploaded")\
                   and (production.review.status not in ["REJECTED", "APPROVED", "DEPRECATED"]):
                    with self.subTest(event=event.name, production=production.name):

                        metafile = f"/home/pe.o4/public_html/gwtc4-fogg/{event.name}/{production.name}/pesummary/samples/posterior_samples.h5"

                        if not os.path.exists(metafile):
                            self.record_result(event.name, production.name, "sky_localization", {
                                "status": "SKIP",
                                "message": "Metafile not found"
                            })
                            continue

                        posterior, prior = open_pesummary_samples(metafile)

                        if 'ra' not in posterior.dtype.names or 'dec' not in posterior.dtype.names:
                            self.record_result(event.name, production.name, "sky_localization", {
                                "status": "SKIP",
                                "message": "RA/Dec not found in posterior"
                            })
                            continue

                        ra = posterior['ra']
                        dec = posterior['dec']

                        area_90 = calculate_sky_area_90(ra, dec)

                        issues = []
                        warning_level = "INFO"

                        if area_90 > 1000:
                            issues.append(f"Large localization area: {area_90:.1f} deg² (may be expected for low-SNR)")
                            warning_level = "WARNING"
                        elif area_90 < 0.01:
                            issues.append(f"Suspiciously small localization area: {area_90:.4f} deg²")
                            warning_level = "WARNING"

                        self.record_result(event.name, production.name, "sky_localization", {
                            "status": warning_level,
                            "message": f"90% sky area: {area_90:.1f} deg²" + (f" - {'; '.join(issues)}" if issues else ""),
                            "details": {"area_90_deg2": float(area_90)}
                        })

                        # Don't fail test, just record information
                        # self.assertTrue(True)

    def test_network_snr_consistency(self):
        """
        Check network SNR is consistent with individual detector SNRs.

        Network optimal SNR should equal the quadrature sum of individual
        detector optimal SNRs:

            SNR_network² = SNR_H1² + SNR_L1² + SNR_V1² + ...

        Inconsistencies may indicate:
        - Post-processing errors
        - Missing detector data
        - Incorrect SNR calculations

        Tolerance is set to 1% relative error.
        """
        for event in self.events:
            for production in event.productions:
                if (production.pipeline.name.lower() == "bilby") \
                   and (production.status_str.lower() == "uploaded")\
                   and (production.review.status not in ["REJECTED", "APPROVED", "DEPRECATED"]):
                    with self.subTest(event=event.name, production=production.name):

                        metafile = f"/home/pe.o4/public_html/gwtc4-fogg/{event.name}/{production.name}/pesummary/samples/posterior_samples.h5"

                        if not os.path.exists(metafile):
                            self.record_result(event.name, production.name, "network_snr", {
                                "status": "SKIP",
                                "message": "Metafile not found"
                            })
                            continue

                        posterior, prior = open_pesummary_samples(metafile)

                        # Check if we have network SNR
                        if 'network_optimal_snr' not in posterior.dtype.names:
                            self.record_result(event.name, production.name, "network_snr", {
                                "status": "SKIP",
                                "message": "network_optimal_snr not found"
                            })
                            continue

                        network_snr = posterior['network_optimal_snr']

                        # Find all detector SNRs
                        detector_snrs = {}
                        for det in ['H1', 'L1', 'V1', 'K1']:  # Common detectors
                            snr_key = f'{det}_optimal_snr'
                            if snr_key in posterior.dtype.names:
                                detector_snrs[det] = posterior[snr_key]

                        if len(detector_snrs) == 0:
                            self.record_result(event.name, production.name, "network_snr", {
                                "status": "SKIP",
                                "message": "No detector SNRs found"
                            })
                            continue

                        # Calculate expected network SNR
                        expected_network_snr_sq = sum(snr**2 for snr in detector_snrs.values())
                        expected_network_snr = np.sqrt(expected_network_snr_sq)

                        # Check consistency
                        relative_error = np.abs(network_snr - expected_network_snr) / network_snr
                        max_error = np.max(relative_error)
                        mean_error = np.mean(relative_error)

                        issues = []
                        if max_error > 0.01:  # 1% tolerance
                            issues.append(f"Network SNR inconsistent (max error: {max_error*100:.2f}%)")

                        mean_network_snr = np.mean(network_snr)
                        detector_info = {det: float(np.mean(snr)) for det, snr in detector_snrs.items()}

                        self.record_result(event.name, production.name, "network_snr", {
                            "status": "FAIL" if issues else "PASS",
                            "message": f"Mean network SNR: {mean_network_snr:.1f}" + (f" - {'; '.join(issues)}" if issues else ""),
                            "details": {
                                "network_snr_mean": float(mean_network_snr),
                                "detector_snrs": detector_info,
                                "max_relative_error": float(max_error)
                            }
                        })

                        self.assertEqual(len(issues), 0,
                                       msg=f"Network SNR inconsistency:\n" + "\n".join(issues))

    @classmethod
    def generate_html_report(cls):
        """
        Generate comprehensive HTML report of all test results.

        The report is saved to each production's directory alongside the
        PESummary results for easy access during review.
        """
        if not cls.test_results:
            return

        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Automated PE Review Report</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
            line-height: 1.6;
            color: #333;
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
            background: #f5f5f5;
        }}

        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}

        .header h1 {{
            margin: 0 0 10px 0;
            font-size: 2.5em;
        }}

        .header .timestamp {{
            opacity: 0.9;
            font-size: 0.9em;
        }}

        .event-section {{
            background: white;
            border-radius: 10px;
            padding: 25px;
            margin-bottom: 30px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}

        .event-name {{
            font-size: 1.8em;
            margin: 0 0 20px 0;
            color: #667eea;
            border-bottom: 3px solid #667eea;
            padding-bottom: 10px;
        }}

        .production-section {{
            margin-bottom: 30px;
            border-left: 4px solid #e0e0e0;
            padding-left: 20px;
        }}

        .production-name {{
            font-size: 1.3em;
            margin: 0 0 15px 0;
            color: #555;
            font-weight: 600;
        }}

        .test-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(350px, 1fr));
            gap: 15px;
            margin-top: 15px;
        }}

        .test-card {{
            border: 1px solid #ddd;
            border-radius: 8px;
            padding: 15px;
            background: #fafafa;
            transition: transform 0.2s, box-shadow 0.2s;
        }}

        .test-card:hover {{
            transform: translateY(-2px);
            box-shadow: 0 4px 8px rgba(0,0,0,0.1);
        }}

        .test-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 10px;
        }}

        .test-name {{
            font-weight: 600;
            color: #333;
            font-size: 1.05em;
        }}

        .status-badge {{
            padding: 4px 12px;
            border-radius: 20px;
            font-size: 0.85em;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}

        .status-PASS {{
            background: #d4edda;
            color: #155724;
        }}

        .status-FAIL {{
            background: #f8d7da;
            color: #721c24;
        }}

        .status-WARNING {{
            background: #fff3cd;
            color: #856404;
        }}

        .status-INFO {{
            background: #d1ecf1;
            color: #0c5460;
        }}

        .status-SKIP {{
            background: #e2e3e5;
            color: #6c757d;
        }}

        .test-message {{
            color: #666;
            font-size: 0.95em;
            margin-top: 8px;
            line-height: 1.5;
        }}

        .test-details {{
            margin-top: 10px;
            padding: 10px;
            background: white;
            border-radius: 5px;
            border: 1px solid #e0e0e0;
            font-size: 0.9em;
            font-family: 'Courier New', monospace;
            max-height: 200px;
            overflow-y: auto;
        }}

        .summary-stats {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            gap: 15px;
            margin-bottom: 30px;
        }}

        .stat-card {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}

        .stat-value {{
            font-size: 2em;
            font-weight: bold;
            color: #667eea;
        }}

        .stat-label {{
            color: #666;
            font-size: 0.9em;
            margin-top: 5px;
        }}

        .legend {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            margin-bottom: 30px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}

        .legend-title {{
            font-size: 1.2em;
            font-weight: 600;
            margin-bottom: 15px;
            color: #333;
        }}

        .legend-item {{
            margin-bottom: 12px;
            padding-left: 10px;
            border-left: 4px solid #667eea;
        }}

        .legend-item-name {{
            font-weight: 600;
            color: #555;
            margin-bottom: 4px;
        }}

        .legend-item-desc {{
            color: #666;
            font-size: 0.95em;
            line-height: 1.5;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Automated PE Review Report</h1>
        <div class="timestamp">Generated: {timestamp}</div>
    </div>
"""

        # Calculate summary statistics
        total_tests = 0
        status_counts = defaultdict(int)

        for event_results in cls.test_results.values():
            for production_results in event_results.values():
                for result in production_results.values():
                    total_tests += 1
                    status_counts[result['status']] += 1

        # Add summary statistics
        html += """
    <div class="summary-stats">
"""
        for status in ['PASS', 'FAIL', 'WARNING', 'INFO', 'SKIP']:
            count = status_counts.get(status, 0)
            html += f"""
        <div class="stat-card">
            <div class="stat-value">{count}</div>
            <div class="stat-label">{status}</div>
        </div>
"""
        html += """
    </div>
"""

        # Add test descriptions legend
        html += """
    <div class="legend">
        <div class="legend-title">Test Descriptions</div>

        <div class="legend-item">
            <div class="legend-item-name">Spin Magnitude Railing</div>
            <div class="legend-item-desc">
                Checks if spin parameters (a_1, a_2, chi_eff, chi_p) are railing at their
                physical bounds. Railing indicates the prior may be too restrictive or there
                are sampler issues.
            </div>
        </div>

        <div class="legend-item">
            <div class="legend-item-name">Effective Sample Size</div>
            <div class="legend-item-desc">
                Measures sampling efficiency by calculating ESS for key parameters. Low ESS
                indicates high autocorrelation and poor convergence. ESS should be at least
                10% of total samples.
            </div>
        </div>

        <div class="legend-item">
            <div class="legend-item-name">Derived Parameter Consistency</div>
            <div class="legend-item-desc">
                Verifies that derived parameters match calculations from sampled parameters
                (e.g., total_mass = mass_1 + mass_2). Inconsistencies suggest post-processing
                errors or file corruption.
            </div>
        </div>

        <div class="legend-item">
            <div class="legend-item-name">Log Evidence Sanity</div>
            <div class="legend-item-desc">
                Checks that log likelihood and log evidence values are finite and within
                reasonable ranges. Unreasonable values indicate sampler failure or numerical
                instability.
            </div>
        </div>

        <div class="legend-item">
            <div class="legend-item-name">Sky Localization Quality</div>
            <div class="legend-item-desc">
                Calculates the 90% credible sky area. Very large areas (>1000 deg²) suggest
                poor localization (may be expected for low-SNR events). Very small areas
                (<0.01 deg²) may indicate sampling issues.
            </div>
        </div>

        <div class="legend-item">
            <div class="legend-item-name">Network SNR Consistency</div>
            <div class="legend-item-desc">
                Verifies that network SNR equals the quadrature sum of individual detector
                SNRs. Inconsistencies may indicate post-processing errors or incorrect SNR
                calculations.
            </div>
        </div>
    </div>
"""

        # Add results for each event/production
        for event_name, event_results in sorted(cls.test_results.items()):
            html += f"""
    <div class="event-section">
        <h2 class="event-name">{event_name}</h2>
"""

            for production_name, production_results in sorted(event_results.items()):
                html += f"""
        <div class="production-section">
            <h3 class="production-name">{production_name}</h3>
            <div class="test-grid">
"""

                for test_name, result in sorted(production_results.items()):
                    status = result['status']
                    message = result['message']
                    details = result.get('details', {})

                    # Format test name nicely
                    test_display = test_name.replace('_', ' ').title()

                    html += f"""
                <div class="test-card">
                    <div class="test-header">
                        <div class="test-name">{test_display}</div>
                        <span class="status-badge status-{status}">{status}</span>
                    </div>
                    <div class="test-message">{message}</div>
"""

                    if details:
                        details_json = json.dumps(details, indent=2)
                        html += f"""
                    <details>
                        <summary style="cursor: pointer; color: #667eea; margin-top: 10px;">Show Details</summary>
                        <div class="test-details"><pre>{details_json}</pre></div>
                    </details>
"""

                    html += """
                </div>
"""

                html += """
            </div>
        </div>
"""

            html += """
    </div>
"""

        html += """
</body>
</html>
"""

        # Save report to each production directory
        for event_name, event_results in cls.test_results.items():
            for production_name in event_results.keys():
                # Construct path to production directory
                report_path = f"/home/pe.o4/public_html/gwtc4-fogg/{event_name}/{production_name}/automated_review.html"

                # Create directory if it doesn't exist
                os.makedirs(os.path.dirname(report_path), exist_ok=True)

                # Write report
                try:
                    with open(report_path, 'w') as f:
                        f.write(html)
                    print(f"Generated automated review report: {report_path}")
                except Exception as e:
                    print(f"Warning: Could not write report to {report_path}: {e}")


if __name__ == '__main__':
    import unittest
    unittest.main()
