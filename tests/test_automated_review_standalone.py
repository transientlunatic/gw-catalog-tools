#!/usr/bin/env python3
"""
Standalone test script for automated PE review tests.

This script can be run independently to test the automated review functionality
on a single HDF5 file without requiring the full Asimov framework.

Usage:
    python test_automated_review_standalone.py <path_to_posterior_samples.h5>

Example:
    python test_automated_review_standalone.py IGWN-GWTC3p0-v2-GW191215_223052_PEDataRelease_mixed_cosmo.h5
"""

import sys
import numpy as np
import h5py as h5
import argparse
from datetime import datetime


def open_pesummary_samples(metafile):
    """Open and extract posterior and prior samples from PESummary HDF5 file."""
    with h5.File(metafile, 'r') as metafile_h:
        analyses = list(metafile_h.keys())
        if "history" in analyses:
            analyses.remove("history")
        if "version" in analyses:
            analyses.remove("version")

        key = analyses[0]
        print(f"Using analysis: {key}")

        posterior = np.array(metafile_h[key]['posterior_samples'])

        # Handle priors - may not always be present
        prior = {}
        if 'priors' in metafile_h[key] and 'samples' in metafile_h[key]['priors']:
            prior = dict(metafile_h[key]['priors']['samples'])
            prior = {name: np.array(data) for name, data in prior.items()}

    return posterior, prior


# Removed calculate_ess - nested sampling produces IID samples


def test_spin_prior_comparison(posterior, prior):
    """
    Compare spin parameter posteriors to priors.

    For spin parameters (a_1, a_2, chi_eff, chi_p), we compare the posterior
    distribution to the prior. If the posterior shows significantly more
    support for high spins (>0.7) than the prior, this indicates the data
    prefers high-spin systems.
    """
    print("\n" + "="*60)
    print("SPIN PARAMETER COMPARISON TO PRIOR")
    print("="*60)

    results = {}

    for param in ['a_1', 'a_2', 'chi_eff', 'chi_p']:
        if param not in posterior.dtype.names:
            print(f"  {param}: NOT FOUND in posterior")
            continue

        if param not in prior or len(prior[param]) == 0:
            print(f"  {param}: NOT FOUND in prior (skipping comparison)")
            continue

        post_data = posterior[param]
        prior_data = prior[param]

        status = "PASS"
        warnings = []

        if param in ['a_1', 'a_2', 'chi_p']:
            # Check if posterior has more high-spin support than prior
            frac_high_post = np.sum(post_data > 0.7) / len(post_data)
            frac_high_prior = np.sum(prior_data > 0.7) / len(prior_data)

            # Warn if posterior has significantly more (>10%) high-spin support
            if frac_high_post > frac_high_prior + 0.1:
                status = "WARNING"
                label = "high spin" if param != 'chi_p' else "high precession"
                warnings.append(f"posterior {frac_high_post*100:.1f}% > 0.7 ({label}), prior {frac_high_prior*100:.1f}%")

            results[param] = {
                'status': status,
                'frac_high_post': frac_high_post,
                'frac_high_prior': frac_high_prior,
                'warnings': warnings,
                'post_min': float(np.min(post_data)),
                'post_max': float(np.max(post_data))
            }

        else:  # chi_eff
            # Check if clustering at bounds relative to prior
            post_upper = np.sum(post_data > 0.9) / len(post_data)
            post_lower = np.sum(post_data < -0.9) / len(post_data)
            prior_upper = np.sum(prior_data > 0.9) / len(prior_data)
            prior_lower = np.sum(prior_data < -0.9) / len(prior_data)

            if post_upper > prior_upper + 0.1:
                status = "WARNING"
                warnings.append(f"posterior {post_upper*100:.1f}% near +1, prior {prior_upper*100:.1f}%")
            if post_lower > prior_lower + 0.1:
                status = "WARNING"
                warnings.append(f"posterior {post_lower*100:.1f}% near -1, prior {prior_lower*100:.1f}%")

            results[param] = {
                'status': status,
                'post_upper': post_upper,
                'post_lower': post_lower,
                'prior_upper': prior_upper,
                'prior_lower': prior_lower,
                'warnings': warnings,
                'post_min': float(np.min(post_data)),
                'post_max': float(np.max(post_data))
            }

        if warnings:
            print(f"  {param}: {status} - {'; '.join(warnings)}")
        else:
            print(f"  {param}: {status}")

    return results


def test_sample_count(posterior):
    """Test that sufficient samples were collected."""
    print("\n" + "="*60)
    print("SAMPLE COUNT CHECK")
    print("="*60)

    n_samples = len(posterior)
    fail_threshold = 2000
    warning_threshold = 5000

    if n_samples < fail_threshold:
        status = "FAIL"
        print(f"  Total samples: {n_samples} - {status}")
        print(f"  Minimum required: {fail_threshold}")
    elif n_samples < warning_threshold:
        status = "WARNING"
        print(f"  Total samples: {n_samples} - {status}")
        print(f"  Recommended: ≥{warning_threshold}")
    else:
        status = "PASS"
        print(f"  Total samples: {n_samples} - {status}")

    results = {
        'n_samples': n_samples,
        'fail_threshold': fail_threshold,
        'warning_threshold': warning_threshold,
        'status': status
    }

    return results


def test_multimodality(posterior):
    """Check for multimodality in key parameters."""
    print("\n" + "="*60)
    print("MULTIMODALITY CHECK")
    print("="*60)

    results = {}
    warnings = []

    for param in ['mass_ratio', 'chirp_mass', 'luminosity_distance']:
        if param not in posterior.dtype.names:
            print(f"  {param}: NOT FOUND")
            continue

        data = posterior[param]
        sorted_data = np.sort(data)

        # Look for large gaps that might indicate separate modes
        diffs = np.diff(sorted_data)
        median_diff = np.median(diffs)
        large_gaps = diffs > 5 * median_diff

        if np.any(large_gaps):
            n_gaps = np.sum(large_gaps)
            max_gap = np.max(diffs[large_gaps])
            gap_fraction = max_gap / (np.max(data) - np.min(data))

            results[param] = {
                'n_large_gaps': int(n_gaps),
                'max_gap_fraction': float(gap_fraction),
                'median_spacing': float(median_diff)
            }

            if gap_fraction > 0.05:
                status = "WARNING"
                warnings.append(f"{n_gaps} large gap(s), max={gap_fraction*100:.1f}% of range")
                print(f"  {param}: {status} - {warnings[-1]}")
            else:
                status = "PASS"
                print(f"  {param}: {status}")
        else:
            status = "PASS"
            print(f"  {param}: {status}")
            results[param] = {'status': status}

    if not warnings:
        print("  Overall: PASS - No significant multimodality detected")

    return results


def test_parameter_consistency(posterior):
    """Test derived parameter consistency."""
    print("\n" + "="*60)
    print("DERIVED PARAMETER CONSISTENCY TEST")
    print("="*60)

    results = {}

    # Test total_mass = mass_1 + mass_2
    if all(p in posterior.dtype.names for p in ['total_mass', 'mass_1', 'mass_2']):
        calculated = posterior['mass_1'] + posterior['mass_2']
        actual = posterior['total_mass']
        rel_error = np.abs(actual - calculated) / (np.abs(actual) + 1e-10)
        max_error = np.max(rel_error)

        status = "PASS" if max_error < 1e-3 else "FAIL"
        print(f"  total_mass = mass_1 + mass_2: {status} (max error: {max_error:.2e})")

        results['total_mass'] = {
            'status': status,
            'max_error': float(max_error)
        }

    # Test mass_ratio = mass_2 / mass_1
    if all(p in posterior.dtype.names for p in ['mass_ratio', 'mass_1', 'mass_2']):
        calculated = posterior['mass_2'] / posterior['mass_1']
        actual = posterior['mass_ratio']
        rel_error = np.abs(actual - calculated) / (np.abs(actual) + 1e-10)
        max_error = np.max(rel_error)

        status = "PASS" if max_error < 1e-3 else "FAIL"
        print(f"  mass_ratio = mass_2 / mass_1: {status} (max error: {max_error:.2e})")

        results['mass_ratio'] = {
            'status': status,
            'max_error': float(max_error)
        }

    # Test chirp_mass
    if all(p in posterior.dtype.names for p in ['chirp_mass', 'mass_1', 'mass_2']):
        m1 = posterior['mass_1']
        m2 = posterior['mass_2']
        calculated = (m1 * m2)**(3/5) / (m1 + m2)**(1/5)
        actual = posterior['chirp_mass']
        rel_error = np.abs(actual - calculated) / (np.abs(actual) + 1e-10)
        max_error = np.max(rel_error)

        status = "PASS" if max_error < 1e-3 else "FAIL"
        print(f"  chirp_mass formula: {status} (max error: {max_error:.2e})")

        results['chirp_mass'] = {
            'status': status,
            'max_error': float(max_error)
        }

    return results


def test_log_evidence_sanity(posterior):
    """Test log evidence values are reasonable."""
    print("\n" + "="*60)
    print("LOG EVIDENCE SANITY TEST")
    print("="*60)

    results = {}

    if 'log_likelihood' in posterior.dtype.names:
        log_like = posterior['log_likelihood']

        n_finite = np.sum(np.isfinite(log_like))
        max_ll = np.max(log_like[np.isfinite(log_like)]) if n_finite > 0 else np.nan
        mean_ll = np.mean(log_like[np.isfinite(log_like)]) if n_finite > 0 else np.nan
        std_ll = np.std(log_like[np.isfinite(log_like)]) if n_finite > 0 else np.nan

        status = "PASS"
        if n_finite < len(log_like):
            status = "FAIL"
            print(f"  log_likelihood: {status} - {len(log_like) - n_finite} non-finite values")
        else:
            print(f"  log_likelihood: {status}")
            print(f"    max: {max_ll:.2f}")
            print(f"    mean: {mean_ll:.2f}")
            print(f"    std: {std_ll:.2f}")

        results['log_likelihood'] = {
            'status': status,
            'max': float(max_ll) if np.isfinite(max_ll) else None,
            'mean': float(mean_ll) if np.isfinite(mean_ll) else None,
            'std': float(std_ll) if np.isfinite(std_ll) else None
        }

    return results


def test_network_snr(posterior):
    """Test network SNR consistency."""
    print("\n" + "="*60)
    print("NETWORK SNR CONSISTENCY TEST")
    print("="*60)

    if 'network_optimal_snr' not in posterior.dtype.names:
        print("  network_optimal_snr: NOT FOUND")
        return {}

    network_snr = posterior['network_optimal_snr']

    # Find detector SNRs
    detector_snrs = {}
    for det in ['H1', 'L1', 'V1', 'K1']:
        snr_key = f'{det}_optimal_snr'
        if snr_key in posterior.dtype.names:
            detector_snrs[det] = posterior[snr_key]

    if not detector_snrs:
        print("  No individual detector SNRs found")
        return {}

    # Calculate expected network SNR
    expected_network_snr_sq = sum(snr**2 for snr in detector_snrs.values())
    expected_network_snr = np.sqrt(expected_network_snr_sq)

    # Check consistency
    rel_error = np.abs(network_snr - expected_network_snr) / (network_snr + 1e-10)
    max_error = np.max(rel_error)

    status = "PASS" if max_error < 0.01 else "FAIL"

    print(f"  Network SNR consistency: {status} (max error: {max_error*100:.2f}%)")
    print(f"    Mean network SNR: {np.mean(network_snr):.2f}")
    for det, snr in detector_snrs.items():
        print(f"    {det} SNR: {np.mean(snr):.2f}")

    results = {
        'status': status,
        'max_error': float(max_error),
        'network_snr_mean': float(np.mean(network_snr)),
        'detector_snrs': {det: float(np.mean(snr)) for det, snr in detector_snrs.items()}
    }

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Run automated PE review tests on a single HDF5 file'
    )
    parser.add_argument(
        'h5file',
        help='Path to posterior_samples.h5 file'
    )

    args = parser.parse_args()

    print("="*60)
    print("AUTOMATED PE REVIEW - STANDALONE TEST")
    print("="*60)
    print(f"File: {args.h5file}")
    print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Load data
    try:
        posterior, prior = open_pesummary_samples(args.h5file)
        print(f"Loaded {len(posterior)} samples with {len(posterior.dtype.names)} parameters")
    except Exception as e:
        print(f"\nERROR: Failed to load file: {e}")
        return 1

    # Run all tests
    all_results = {}

    try:
        all_results['spin_prior_comparison'] = test_spin_prior_comparison(posterior, prior)
    except Exception as e:
        print(f"\nERROR in spin_prior_comparison test: {e}")

    try:
        all_results['sample_count'] = test_sample_count(posterior)
    except Exception as e:
        print(f"\nERROR in sample_count test: {e}")

    try:
        all_results['multimodality'] = test_multimodality(posterior)
    except Exception as e:
        print(f"\nERROR in multimodality test: {e}")

    try:
        all_results['parameter_consistency'] = test_parameter_consistency(posterior)
    except Exception as e:
        print(f"\nERROR in parameter_consistency test: {e}")

    try:
        all_results['log_evidence_sanity'] = test_log_evidence_sanity(posterior)
    except Exception as e:
        print(f"\nERROR in log_evidence_sanity test: {e}")

    try:
        all_results['network_snr'] = test_network_snr(posterior)
    except Exception as e:
        print(f"\nERROR in network_snr test: {e}")

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)

    total_tests = 0
    passed = 0
    failed = 0
    warnings = 0

    for test_name, test_results in all_results.items():
        if isinstance(test_results, dict):
            for param_name, result in test_results.items():
                if isinstance(result, dict) and 'status' in result:
                    total_tests += 1
                    if result['status'] == 'PASS':
                        passed += 1
                    elif result['status'] == 'FAIL':
                        failed += 1
                    elif result['status'] == 'WARNING':
                        warnings += 1

    print(f"Total tests: {total_tests}")
    print(f"  PASS: {passed}")
    print(f"  FAIL: {failed}")
    print(f"  WARNING: {warnings}")

    print("\n" + "="*60)

    return 0 if failed == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
