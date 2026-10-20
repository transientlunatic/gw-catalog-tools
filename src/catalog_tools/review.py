"""
The checks of the automated PE review.

Everything which is common to the asimov-driven tests and the command line
lives here so that the two cannot drift apart.  This module needs only numpy,
scipy and h5py, and depends on asimov only to locate files through a project's
config, which is imported lazily.

Environment variables
---------------------
PE_REVIEW_INCLUDE_REVIEWED
    Set to ``1`` to run the tests on analyses which have already been reviewed
    (by default only unreviewed analyses are tested).  Useful for validating
    the tests against past decisions.
PE_REVIEW_REPORT_DIR
    Write the HTML reports under this directory (as
    ``<dir>/<event>/<analysis>/automated_review.html``) instead of alongside
    the PESummary pages in the web root.
"""

import os

import h5py as h5
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

REVIEWED_STATES = ("REJECTED", "APPROVED", "DEPRECATED")

# Network SNR columns, in order of preference
SNR_KEYS = ("network_matched_filter_snr", "network_optimal_snr")

# Thresholds for the plausibility checks.  These are tuned for a catalogue of
# near-threshold candidates, where a genuine signal has a network SNR of
# roughly 6--15.
SNR_NOISE_THRESHOLD = 4.0
SNR_LOUD_THRESHOLD = 20.0
LN_BAYES_FACTOR_THRESHOLD = 1.0
DETECTOR_DOMINANCE_THRESHOLD = 0.75
NEFF_FAIL_THRESHOLD = 1000
NEFF_WARNING_THRESHOLD = 2500


# ---------------------------------------------------------------------------
# Selecting analyses
# ---------------------------------------------------------------------------

def include_reviewed():
    """True if reviewed analyses should be tested as well."""
    return os.environ.get("PE_REVIEW_INCLUDE_REVIEWED", "").lower() in ("1", "true", "yes")


def is_bilby(production):
    """True for a bilby analysis.

    The pipeline name is "Bilby" in asimov 0.6 and "bilby" in 0.7 onwards, so
    it is always compared case-insensitively.
    """
    return production.pipeline.name.lower() == "bilby"


def is_uploaded(production):
    return production.status_str.lower() == "uploaded"


def is_approved(production):
    status = production.review.status
    return status is not None and status.upper() == "APPROVED"


def should_review(production):
    """True for uploaded bilby analyses which still need to be reviewed.

    Analyses which already have a review are skipped unless
    ``PE_REVIEW_INCLUDE_REVIEWED`` is set.
    """
    if not (is_bilby(production) and is_uploaded(production)):
        return False
    return include_reviewed() or production.review.status not in REVIEWED_STATES


# ---------------------------------------------------------------------------
# Locating and reading results
# ---------------------------------------------------------------------------

def metafile_candidates(event_name, production_name, webroot, project_root, storage):
    """
    Where an analysis' PESummary samples file might be, in order of preference:

    1. in the web root, where PESummary writes it;
    2. in the results store of the project in the current directory;
    3. in the results store of the project which the config names.

    The second is tried before the third because a project's config can name a
    different copy of the project (an asimov project can be moved or copied
    while its config still holds the old absolute paths), and the analyses of
    the project being tested are in the current directory.
    """
    return [
        # os.path.join keeps an absolute webroot as it is
        os.path.join(project_root, webroot, event_name, production_name,
                     "pesummary", "samples", "posterior_samples.h5"),
        os.path.join(storage, event_name, production_name, "posterior_samples.h5"),
        os.path.join(project_root, storage, event_name, production_name,
                     "posterior_samples.h5"),
    ]


def candidate_metafiles(event_name, production_name):
    """:func:`metafile_candidates`, with the locations taken from the project config."""
    from asimov import config

    return metafile_candidates(
        event_name, production_name,
        webroot=config.get("general", "webroot"),
        project_root=config.get("project", "root", fallback=""),
        storage=config.get("storage", "directory", fallback="results"),
    )


def get_metafile_path(event_name, production_name):
    """Path to the posterior samples file for an analysis.

    Returns the first of :func:`metafile_candidates` which exists.  If none
    does the first is returned, so that the caller can report it.
    """
    candidates = candidate_metafiles(event_name, production_name)
    for path in candidates:
        if os.path.exists(path):
            return path
    return candidates[0]


def get_report_path(event_name, production_name):
    """Where to write the HTML report for an analysis."""
    root = os.environ.get("PE_REVIEW_REPORT_DIR")
    if root:
        return os.path.join(root, event_name, production_name, "automated_review.html")

    from asimov import config

    webroot = config.get("general", "webroot")
    project_root = config.get("project", "root", fallback="")
    return os.path.join(project_root, webroot, event_name, production_name,
                        "automated_review.html")


def _first_analysis(metafile_h):
    keys = [key for key in metafile_h.keys() if key not in ("history", "version")]
    return keys[0]


def open_pesummary_samples(metafile):
    """Posterior samples (structured array) and prior samples (dict of arrays)."""
    with h5.File(metafile, "r") as metafile_h:
        key = _first_analysis(metafile_h)
        posterior = np.array(metafile_h[key]["posterior_samples"])
        prior = {}
        if "priors" in metafile_h[key] and "samples" in metafile_h[key]["priors"]:
            prior = {name: np.array(data)
                     for name, data in metafile_h[key]["priors"]["samples"].items()}
    return posterior, prior


def read_run_info(metafile):
    """Sampler statistics stored in the file; entries which are missing are None."""
    info = {"ln_bayes_factor": None, "ln_evidence": None, "neffsamples": None}
    with h5.File(metafile, "r") as metafile_h:
        meta = metafile_h[_first_analysis(metafile_h)].get("meta_data", {})
        sampler = meta.get("sampler", {})
        for name in ("ln_bayes_factor", "ln_evidence"):
            if name in sampler:
                info[name] = float(np.array(sampler[name]).ravel()[0])
        stats = meta.get("other", {}).get("run_statistics", {})
        if "neffsamples" in stats:
            info["neffsamples"] = float(np.array(stats["neffsamples"]).ravel()[0])
    return info


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

def find_modes(samples, log=False, min_prominence=0.08, min_mass=0.03,
               max_samples=20000, seed=0):
    """
    Find the modes of a one-dimensional posterior.

    The samples are smoothed with a Gaussian kernel (Scott's bandwidth), and a
    mode is a peak of the smoothed density which

    * rises above its surroundings by at least ``min_prominence`` of the
      highest peak (so that ripples and tails are ignored), and
    * contains at least ``min_mass`` of the posterior probability, where
      the mass is shared between peaks at the density minimum between them.

    The default prominence was chosen on synthetic samples: no false positives
    in 100 draws each from Gaussian, log-normal, gamma and bounded
    distributions (1 in 100 for a perfectly flat one), and reliable detection
    of a second mode carrying 20% or more of the probability.

    This replaces an earlier test which looked for gaps in the sorted samples.
    That compared gaps to the median spacing in the dense core, so the
    sparsely-sampled tails of any posterior looked like gaps, and it could not
    see two modes which overlap and so have no gap at all.

    Parameters
    ----------
    samples : array_like
    log : bool
        Estimate the density in log(samples), for strictly positive,
        scale-like parameters such as distance.

    Returns
    -------
    list of dict
        One ``{"location", "mass", "height"}`` entry per mode, ordered by
        location; ``height`` is relative to the highest peak.
    """
    x = np.asarray(samples, dtype=float)
    x = x[np.isfinite(x)]
    if log:
        x = np.log(x[x > 0])
    if x.size < 100 or np.std(x) == 0:
        return []
    if x.size > max_samples:
        x = np.random.default_rng(seed).choice(x, max_samples, replace=False)

    bandwidth = x.std() * x.size ** (-1 / 5)
    low, high = x.min() - 3 * bandwidth, x.max() + 3 * bandwidth
    nbins = int(np.clip((high - low) / bandwidth * 10, 256, 4000))
    counts, edges = np.histogram(x, bins=nbins, range=(low, high))
    width = edges[1] - edges[0]
    density = gaussian_filter1d(counts.astype(float), bandwidth / width, mode="constant")
    centres = 0.5 * (edges[1:] + edges[:-1])

    # Padding with zeros lets a peak against the edge of the samples (such as
    # at a prior bound) count as a peak.
    peaks, _ = find_peaks(np.pad(density, 1), prominence=min_prominence * density.max())
    peaks = peaks - 1
    if len(peaks) == 0:
        return []

    cuts = [0]
    for left, right in zip(peaks[:-1], peaks[1:]):
        cuts.append(left + int(np.argmin(density[left:right + 1])))
    cuts.append(len(density))

    total = density.sum()
    modes = []
    for index, peak in enumerate(peaks):
        mass = density[cuts[index]:cuts[index + 1]].sum() / total
        if mass >= min_mass:
            location = centres[peak]
            modes.append({
                "location": float(np.exp(location) if log else location),
                "mass": float(mass),
                "height": float(density[peak] / density.max()),
            })
    return modes


def check_parameter_consistency(posterior, param_derived, param_inputs, calculation_func,
                                tolerance=1e-3):
    """
    Check that a derived parameter matches its calculation from other parameters.

    Returns ``(is_consistent, max_relative_error)``, or ``(None, None)`` if any
    of the parameters is missing from the posterior.
    """
    for param in [param_derived] + param_inputs:
        if param not in posterior.dtype.names:
            return None, None

    derived_values = posterior[param_derived]
    calculated = calculation_func(*[posterior[param] for param in param_inputs])
    relative_errors = np.abs(derived_values - calculated) / (np.abs(derived_values) + 1e-10)
    max_error = np.max(relative_errors)
    return bool(max_error < tolerance), float(max_error)


def calculate_sky_area_90(ra, dec, n_ra=64, n_dec=32):
    """
    90% credible area of the sky localisation, in square degrees.

    The samples are binned on a grid which is uniform in right ascension and
    in sin(declination), so every cell has the same solid angle, and the area
    is the number of cells needed to enclose 90% of the samples.  The result
    is limited by the cell size (about 20 deg^2 with the default grid).

    Returns
    -------
    area : float
        The 90% area in square degrees.
    cell_area : float
        The area of one cell, in square degrees.
    """
    ra = np.mod(np.asarray(ra, dtype=float), 2 * np.pi)
    sin_dec = np.sin(np.asarray(dec, dtype=float))
    counts, _, _ = np.histogram2d(ra, sin_dec, bins=[n_ra, n_dec],
                                  range=[[0, 2 * np.pi], [-1, 1]])
    cell_area = (2 * np.pi / n_ra) * (2 / n_dec) * (180 / np.pi) ** 2
    cumulative = np.cumsum(np.sort(counts.ravel())[::-1]) / counts.sum()
    n_cells = int(np.searchsorted(cumulative, 0.90)) + 1
    return n_cells * cell_area, cell_area


def _sigmoid(support, steepness=0.03):
    """Empirical sigmoid used to size the suggested new bound."""
    return 0.66 / (1 + np.exp(-steepness * (support - 30)))


def railing_check(samples, nbins=50, tolerance=2.0, label=None, input_min_edge=None):
    """
    Check whether a posterior is railing against the prior.

    The posterior is histogrammed, and it is railing at an end if the bin at
    that end holds more than ``tolerance`` percent of the density of the
    highest bin.  The criterion, and the sizing of the suggested new bound, are
    those of the railing check in peconfigurator (itself taken from LAPLACE),
    which this reproduces without needing peconfigurator and its dependencies;
    the tests compare the two when peconfigurator is installed.  The one
    difference is that the suggested upper bound is sized from the support at
    the upper end, where peconfigurator uses the support at the lower end.

    Parameters
    ----------
    samples : array_like
    nbins : int
    tolerance : float
        Percent of the highest bin above which an end bin counts as railing.
    label : str, optional
        "mass_ratio" is never reported as railing at the upper end, since
        q = 1 is a genuine boundary.
    input_min_edge : float, optional
        The lower edge of the prior.  Railing at the lower end is ignored if
        the prior already extends below the suggested new bound.  This avoids
        false positives for distance, which can legitimately pile up near zero.

    Returns
    -------
    railing_low, railing_high : bool
    recommendation : str
    suggested_bounds : tuple
        The suggested ``(lower, upper)`` bounds; ``None`` where there is no
        railing.
    """
    samples = np.asarray(samples, dtype=float)
    hist, _ = np.histogram(samples, bins=nbins, density=True)
    highest = np.amax(hist)
    lower_support = hist[0] / highest * 100
    higher_support = hist[-1] / highest * 100

    railing_low = bool(lower_support > tolerance)
    railing_high = bool(higher_support > tolerance)
    if label == "mass_ratio":
        railing_high = False

    bound_low = bound_high = None
    recommendation = ""
    if railing_low:
        edge = np.round(np.amin(samples), 4)
        bound_low = np.round((1 - _sigmoid(lower_support)) * edge, 4)
        recommendation += f"Railing at lower end at {edge}. Suggested new bound is {bound_low}. "
        if input_min_edge is not None and input_min_edge < bound_low:
            railing_low, bound_low, recommendation = False, None, ""
    if railing_high:
        edge = np.round(np.amax(samples), 4)
        bound_high = np.round((1 + _sigmoid(higher_support)) * edge, 4)
        recommendation += f"Railing at higher end at {edge}. Suggested new bound is {bound_high}"

    return railing_low, railing_high, recommendation, (bound_low, bound_high)


# ---------------------------------------------------------------------------
# The checks
#
# Each takes the data it needs and returns a result dictionary with a
# ``status`` (PASS, WARNING, FAIL, INFO or SKIP), a one-line ``message`` and a
# ``details`` dictionary.  The asimov tests and the standalone script both
# call these, so that they always agree.
# ---------------------------------------------------------------------------

def _result(status, message, **details):
    return {"status": status, "message": message, "details": details}


def check_spin_prior(posterior, prior):
    """
    Compare spin posteriors to the prior.

    For a_1, a_2 and chi_p the fraction of samples above 0.7 is compared with
    that of the prior, and for chi_eff the fractions beyond +-0.9.  A posterior
    with more than 10% more support there than the prior is pinned against the
    spin bounds, which is characteristic of a glitch or of a noise-dominated
    analysis.  Screening only: this is a WARNING, never a FAIL.
    """
    warnings = []
    for param in ("a_1", "a_2", "chi_eff", "chi_p"):
        if param not in posterior.dtype.names or len(prior.get(param, [])) == 0:
            continue
        if param == "chi_eff":
            for label, test in (("+1", lambda x: x > 0.9), ("-1", lambda x: x < -0.9)):
                post_frac = np.mean(test(posterior[param]))
                prior_frac = np.mean(test(prior[param]))
                if post_frac > prior_frac + 0.1:
                    warnings.append(f"chi_eff: {post_frac*100:.1f}% near {label}, "
                                    f"prior: {prior_frac*100:.1f}%")
        else:
            post_frac = np.mean(posterior[param] > 0.7)
            prior_frac = np.mean(prior[param] > 0.7)
            if post_frac > prior_frac + 0.1:
                label = "high precession" if param == "chi_p" else "high spin"
                warnings.append(f"{param}: {post_frac*100:.1f}% > 0.7 ({label}), "
                                f"prior: {prior_frac*100:.1f}%")
    if warnings:
        return _result("WARNING", "; ".join(warnings), warnings=warnings)
    return _result("PASS", "Spins are not pushed towards the bounds relative to the prior")


def check_effective_sample_size(n_samples, info):
    """
    Check the effective number of posterior samples is sufficient.

    The sampler's own estimate (``run_statistics/neffsamples``) is used when it
    is in the file, otherwise the number of samples, which is an upper bound.
    FAIL below NEFF_FAIL_THRESHOLD, WARNING below NEFF_WARNING_THRESHOLD.
    """
    neff = info.get("neffsamples")
    source = "sampler estimate"
    if neff is None:
        neff, source = float(n_samples), "sample count"
    if neff < NEFF_FAIL_THRESHOLD:
        status = "FAIL"
    elif neff < NEFF_WARNING_THRESHOLD:
        status = "WARNING"
    else:
        status = "PASS"
    return _result(status, f"{neff:.0f} effective samples ({source}) from {n_samples} samples",
                   neff=neff, n_samples=n_samples, source=source,
                   fail_threshold=NEFF_FAIL_THRESHOLD,
                   warning_threshold=NEFF_WARNING_THRESHOLD)


def check_multimodality(posterior):
    """
    Flag parameters whose posterior has more than one significant mode.

    See :func:`find_modes`.  Luminosity distance is examined in log space.
    theta_JN is not checked, since its bimodality about pi/2 is expected for a
    two-detector network.  Screening only: a WARNING, never a FAIL.
    """
    warnings = []
    details = {}
    for param, log in (("mass_ratio", False), ("chirp_mass", False),
                       ("luminosity_distance", True)):
        if param not in posterior.dtype.names:
            continue
        modes = find_modes(posterior[param], log=log)
        details[param] = modes
        if len(modes) > 1:
            described = ", ".join(f"{m['location']:.3g} ({m['mass']*100:.0f}%)" for m in modes)
            warnings.append(f"{param}: {len(modes)} modes at {described}")
    if warnings:
        return _result("WARNING", "Multiple modes: " + "; ".join(warnings), **details)
    return _result("PASS", "Each parameter has a single mode", **details)


def check_derived_parameters(posterior):
    """Check total mass, mass ratio and chirp mass against the component masses (0.1%)."""
    def chirp_mass(m1, m2):
        return (m1 * m2) ** (3 / 5) / (m1 + m2) ** (1 / 5)

    relations = (
        ("total_mass", ["mass_1", "mass_2"], lambda m1, m2: m1 + m2, "mass_1 + mass_2"),
        ("mass_ratio", ["mass_2", "mass_1"], lambda m2, m1: m2 / m1, "mass_2 / mass_1"),
        ("chirp_mass", ["mass_1", "mass_2"], chirp_mass, "the component masses"),
    )
    issues, checks = [], []
    for derived, inputs, function, description in relations:
        consistent, max_error = check_parameter_consistency(posterior, derived, inputs, function)
        if consistent is None:
            continue
        checks.append(f"{derived} consistency: {max_error:.2e}")
        if not consistent:
            issues.append(f"{derived} inconsistent with {description} (max error: {max_error:.2e})")
    if issues:
        return _result("FAIL", "; ".join(issues), checks=checks, issues=issues)
    return _result("PASS", "All checked parameters consistent", checks=checks, issues=issues)


def check_likelihood(posterior):
    """
    Check the log-likelihood and log-prior samples are finite, and that the
    log-likelihood has some spread (a nearly constant likelihood means the
    likelihood surface was not explored).
    """
    issues, details = [], {}
    if "log_likelihood" in posterior.dtype.names:
        log_like = posterior["log_likelihood"]
        if not np.all(np.isfinite(log_like)):
            issues.append(f"log_likelihood has {np.sum(~np.isfinite(log_like))} non-finite values")
        else:
            details["log_likelihood"] = {"max": float(np.max(log_like)),
                                         "mean": float(np.mean(log_like)),
                                         "std": float(np.std(log_like))}
            if np.std(log_like) < 0.1:
                issues.append(f"log_likelihood has very small std: {np.std(log_like):.2e}")
    if "log_prior" in posterior.dtype.names:
        log_prior = posterior["log_prior"]
        if not np.all(np.isfinite(log_prior)):
            issues.append(f"log_prior has {np.sum(~np.isfinite(log_prior))} non-finite values")
    if issues:
        return _result("FAIL", "; ".join(issues), **details)
    return _result("PASS", "Log-likelihood and log-prior values are reasonable", **details)


def check_snr_plausibility(posterior, info):
    """
    Check the signal strength is plausible for a near-threshold candidate.

    * ln Bayes factor below LN_BAYES_FACTOR_THRESHOLD, or median network
      matched-filter SNR below SNR_NOISE_THRESHOLD: the data do not support a
      signal and the posterior is dominated by the prior.
    * Median network SNR above SNR_LOUD_THRESHOLD: far louder than a candidate
      below the significance threshold should be, which usually means the
      analysis has locked on to a loud noise transient.
    * More than DETECTOR_DOMINANCE_THRESHOLD of the SNR^2 in one detector (of
      two or more): the candidate may be a glitch in that detector.

    The details give each of the three conditions separately, as
    ``consistent_with_noise``, ``implausibly_loud`` and
    ``single_detector_dominated``.  Screening only: a WARNING, never a FAIL.
    """
    names = posterior.dtype.names
    snr_key = next((k for k in SNR_KEYS if k in names), None)
    snr = float(np.median(np.real(posterior[snr_key]))) if snr_key else None
    ln_bf = info.get("ln_bayes_factor")
    warnings = []

    noise_like = bool((ln_bf is not None and ln_bf < LN_BAYES_FACTOR_THRESHOLD)
                      or (snr is not None and snr < SNR_NOISE_THRESHOLD))
    if ln_bf is not None and ln_bf < LN_BAYES_FACTOR_THRESHOLD:
        warnings.append(f"ln Bayes factor {ln_bf:.2f} < {LN_BAYES_FACTOR_THRESHOLD:g}: "
                        "consistent with noise")
    if snr is not None and snr < SNR_NOISE_THRESHOLD:
        warnings.append(f"median network SNR {snr:.1f} < {SNR_NOISE_THRESHOLD:g}: "
                        "consistent with noise")

    loud = bool(snr is not None and snr > SNR_LOUD_THRESHOLD)
    if loud:
        warnings.append(f"median network SNR {snr:.1f} > {SNR_LOUD_THRESHOLD:g}: "
                        "implausibly loud for a near-threshold candidate")

    details = {"median_network_snr": snr, "ln_bayes_factor": ln_bf,
               "consistent_with_noise": noise_like, "implausibly_loud": loud,
               "single_detector_dominated": False}
    detector_snr_sq = {det: float(np.median(posterior[f"{det}_optimal_snr"])) ** 2
                       for det in ("H1", "L1", "V1", "K1") if f"{det}_optimal_snr" in names}
    if len(detector_snr_sq) > 1:
        dominant = max(detector_snr_sq, key=detector_snr_sq.get)
        fraction = detector_snr_sq[dominant] / sum(detector_snr_sq.values())
        details["dominant_detector"] = {"detector": dominant, "snr_squared_fraction": fraction}
        if fraction > DETECTOR_DOMINANCE_THRESHOLD:
            details["single_detector_dominated"] = True
            warnings.append(f"{dominant} contributes {fraction*100:.0f}% of the SNR^2")

    if warnings:
        return _result("WARNING", "; ".join(warnings), **details)
    summary = f"Network SNR {snr:.1f}" if snr is not None else "Network SNR unavailable"
    if ln_bf is not None:
        summary += f", ln BF {ln_bf:.1f}"
    return _result("PASS", summary, **details)


def check_sky_localization(posterior):
    """
    Record the 90% credible sky area.

    Areas of thousands of square degrees are normal for near-threshold
    candidates, so this is informational.  A WARNING is raised only if the whole
    90% region fits in a couple of grid cells, which suggests the sky samples
    have collapsed.
    """
    if "ra" not in posterior.dtype.names or "dec" not in posterior.dtype.names:
        return _result("SKIP", "RA/Dec not found in posterior")
    area, cell_area = calculate_sky_area_90(posterior["ra"], posterior["dec"])
    collapsed = area <= 2 * cell_area
    message = f"90% sky area: {area:.0f} deg²"
    if collapsed:
        message += " - sky samples are concentrated in almost a single cell"
    return _result("WARNING" if collapsed else "INFO", message,
                   area_90_deg2=float(area), cell_area_deg2=float(cell_area))


def check_network_snr(posterior):
    """The network optimal SNR should be the quadrature sum of the detector SNRs (1%)."""
    names = posterior.dtype.names
    if "network_optimal_snr" not in names:
        return _result("SKIP", "network_optimal_snr not found")
    detector_snrs = {det: posterior[f"{det}_optimal_snr"]
                     for det in ("H1", "L1", "V1", "K1") if f"{det}_optimal_snr" in names}
    if not detector_snrs:
        return _result("SKIP", "No detector SNRs found")

    network = posterior["network_optimal_snr"]
    expected = np.sqrt(sum(snr ** 2 for snr in detector_snrs.values()))
    max_error = float(np.max(np.abs(network - expected) / np.maximum(network, 1e-10)))
    details = {"network_snr_mean": float(np.mean(network)),
               "detector_snrs": {det: float(np.mean(snr)) for det, snr in detector_snrs.items()},
               "max_relative_error": max_error}
    message = f"Mean network SNR: {np.mean(network):.1f}"
    if max_error > 0.01:
        return _result("FAIL", message + f" - Network SNR inconsistent (max error: {max_error*100:.2f}%)",
                       **details)
    return _result("PASS", message, **details)


def check_railing(posterior, prior, param):
    """
    Check one parameter's posterior for railing against the prior bounds.

    Used for chirp mass and luminosity distance.  Railing means the prior is
    cutting off real support and should be widened, so this is a FAIL, and the
    result includes the suggested new bounds.  For distance the prior minimum
    is passed on so that piling up towards zero is not reported.
    """
    if param not in posterior.dtype.names:
        return _result("SKIP", f"{param} not found")
    kwargs = {}
    if param == "luminosity_distance" and len(prior.get(param, [])):
        kwargs["input_min_edge"] = float(np.min(prior[param]))
    low, high, _, (bound_low, bound_high) = railing_check(posterior[param], label=param, **kwargs)
    messages = []
    if low:
        messages.append(f"railing at the lower bound; recommend {bound_low}")
    if high:
        messages.append(f"railing at the upper bound; recommend {bound_high}")
    details = {"railing_lower": low, "railing_upper": high,
               "suggested_lower_bound": None if bound_low is None else float(bound_low),
               "suggested_upper_bound": None if bound_high is None else float(bound_high)}
    if messages:
        return _result("FAIL", f"{param}: " + "; ".join(messages), **details)
    return _result("PASS", f"No railing in {param}", **details)


def run_all_checks(posterior, prior, info):
    """Every check, as ``{name: result}``."""
    return {
        "chirp_mass_railing": check_railing(posterior, prior, "chirp_mass"),
        "distance_railing": check_railing(posterior, prior, "luminosity_distance"),
        "spin_prior_comparison": check_spin_prior(posterior, prior),
        "effective_sample_size": check_effective_sample_size(len(posterior), info),
        "multimodality": check_multimodality(posterior),
        "parameter_consistency": check_derived_parameters(posterior),
        "likelihood_sanity": check_likelihood(posterior),
        "snr_plausibility": check_snr_plausibility(posterior, info),
        "sky_localization": check_sky_localization(posterior),
        "network_snr": check_network_snr(posterior),
    }
