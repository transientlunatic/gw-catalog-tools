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


def calculate_ess(samples):
    """Calculate effective sample size using autocorrelation."""
    n = len(samples)
    if n < 2:
        return n

    samples = samples - np.mean(samples)
    c0 = np.dot(samples, samples) / n
    if c0 == 0:
        return n

    autocorr_sum = 0
    for k in range(1, min(n, 1000)):
        ck = np.dot(samples[:-k], samples[k:]) / (n - k)
        autocorr_sum += ck / c0
        if ck / c0 < 0.05:
            break

    ess = n / (1 + 2 * autocorr_sum)
    return max(1, ess)


def test_spin_railing(posterior):
    """
    Test for spin parameter railing (simplified without peconfigurator).

    Railing at 0 or 1 is acceptable for spin magnitudes (a_1, a_2, chi_p).
    Warns if >10% of samples are above 0.7 (high-spin systems).
    For chi_eff, check both bounds as it can be negative.
    """
    print("\n" + "="*60)
    print("SPIN MAGNITUDE RAILING TEST")
    print("="*60)

    results = {}

    for param in ['a_1', 'a_2', 'chi_eff', 'chi_p']:
        if param not in posterior.dtype.names:
            print(f"  {param}: NOT FOUND")
            continue

        data = posterior[param]

        status = "PASS"
        warnings = []

        if param in ['a_1', 'a_2', 'chi_p']:
            # For spin magnitudes: don't flag railing at 0 or 1 (both are physical)
            # Instead, warn if significant support above 0.7
            frac_high_spin = np.sum(data > 0.7) / len(data)

            if frac_high_spin > 0.1:
                status = "WARNING"
                label = "high spin" if param != 'chi_p' else "high precession"
                warnings.append(f"{frac_high_spin*100:.1f}% of samples > 0.7 ({label})")

            results[param] = {
                'status': status,
                'frac_high_spin': frac_high_spin,
                'warnings': warnings,
                'min': float(np.min(data)),
                'max': float(np.max(data))
            }

        else:  # chi_eff
            # chi_eff can be negative, check railing at both bounds
            data_range = 2.0  # -1 to 1
            threshold = 0.05 * data_range

            n_near_lower = np.sum(data < (-1.0 + threshold))
            n_near_upper = np.sum(data > (1.0 - threshold))

            frac_lower = n_near_lower / len(data)
            frac_upper = n_near_upper / len(data)

            issues = []
            if frac_upper > 0.1:
                status = "WARNING"
                issues.append(f"upper bound ({frac_upper*100:.1f}%)")

            if frac_lower > 0.1:
                status = "WARNING"
                issues.append(f"lower bound ({frac_lower*100:.1f}%)")

            results[param] = {
                'status': status,
                'frac_lower': frac_lower,
                'frac_upper': frac_upper,
                'issues': issues,
                'min': float(np.min(data)),
                'max': float(np.max(data))
            }

            warnings = issues

        if warnings:
            print(f"  {param}: {status} - {'; '.join(warnings)}")
        else:
            print(f"  {param}: {status}")

    return results


def test_effective_sample_size(posterior):
    """Test effective sample size for key parameters."""
    print("\n" + "="*60)
    print("EFFECTIVE SAMPLE SIZE TEST")
    print("="*60)

    n_samples = len(posterior)
    print(f"  Total samples: {n_samples}")

    results = {}
    params_to_check = ['luminosity_distance', 'chirp_mass', 'mass_ratio', 'mass_1', 'mass_2']

    for param in params_to_check:
        if param not in posterior.dtype.names:
            print(f"  {param}: NOT FOUND")
            continue

        ess = calculate_ess(posterior[param])
        ess_fraction = ess / n_samples

        status = "PASS" if ess_fraction >= 0.1 else "WARNING"

        print(f"  {param}: ESS={ess:.0f} ({ess_fraction*100:.1f}% of total) - {status}")

        results[param] = {
            'ess': float(ess),
            'fraction': float(ess_fraction),
            'status': status
        }

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
        all_results['spin_railing'] = test_spin_railing(posterior)
    except Exception as e:
        print(f"\nERROR in spin_railing test: {e}")

    try:
        all_results['effective_sample_size'] = test_effective_sample_size(posterior)
    except Exception as e:
        print(f"\nERROR in effective_sample_size test: {e}")

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
