"""
Automated PE Review Tests
==========================

This module provides automated quality checks for parameter estimation
analyses. Tests run once per production and generate an HTML report for review.

The tests check for common failure modes including:
- Posterior support for high spins compared with the prior
- Poor sampling quality (low effective sample size)
- Multiple modes in the key mass and distance parameters
- Derived parameter inconsistencies
- Likelihood sanity
- Implausible signal strength (consistent with noise, or too loud for a
  near-threshold candidate) and strong dominance by a single detector
- Sky localization quality
- Network SNR consistency

The helpers are in ``review_utils``, which is shared with the standalone
script.  By default only uploaded bilby analyses which have not yet been
reviewed are tested; set ``PE_REVIEW_INCLUDE_REVIEWED=1`` to test the others,
and ``PE_REVIEW_REPORT_DIR`` to write the reports somewhere other than the
web root.  Tests are designed to run without requiring trained ML models.
"""

from asimov.testing import AsimovTest
import os
from datetime import datetime
from collections import defaultdict
import json

import review_utils as utils
from review_utils import (
    get_metafile_path,
    get_report_path,
    open_pesummary_samples,
    read_run_info,
    should_review,
)


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

    def analyses_to_review(self, test_key):
        """
        The (event, production, metafile) triples to be tested.

        A SKIP is recorded for any analysis whose samples can't be found, and
        the calling test is skipped if there is nothing to test at all, so
        that a test which has checked nothing can't be mistaken for a pass.
        """
        todo = []
        for event in self.events:
            for production in event.productions:
                if not should_review(production):
                    continue
                metafile = get_metafile_path(event.name, production.name)
                if not os.path.exists(metafile):
                    self.record_result(event.name, production.name, test_key, {
                        "status": "SKIP",
                        "message": f"Samples not found at {metafile}"
                    })
                    continue
                todo.append((event, production, metafile))
        if not todo:
            self.skipTest("No analyses to review")
        return todo

    def run_check(self, test_key, check, blocking=False):
        """
        Run ``check`` on every analysis to be reviewed and record the results.

        ``check`` is called with the keyword arguments ``posterior``, ``prior``
        and ``info`` (the sampler statistics) and returns a result dictionary.
        If ``blocking`` is set a FAIL fails the test; otherwise results are only
        recorded in the report, as for screening checks.

        Failures are collected and reported together, each with the event and
        analysis it belongs to, so that a failure is never ambiguous.
        """
        failures = []
        for event, production, metafile in self.analyses_to_review(test_key):
            label = f"{event.name} {production.name}"
            try:
                posterior, prior = open_pesummary_samples(metafile)
                result = check(posterior=posterior, prior=prior, info=read_run_info(metafile))
            except Exception as error:
                self.record_result(event.name, production.name, test_key, {
                    "status": "FAIL", "message": f"Could not run check: {error!r}"})
                failures.append(f"{label}: could not run check: {error!r}")
                continue
            self.record_result(event.name, production.name, test_key, result)
            if blocking and result["status"] == "FAIL":
                failures.append(f"{label}: {result['message']}")
        if failures:
            self.fail("\n".join(failures))

    def test_spin_prior_comparison(self):
        """Spins should not be pushed towards their bounds relative to the prior (warning only)."""
        self.run_check("spin_prior_comparison",
                       lambda posterior, prior, info: utils.check_spin_prior(posterior, prior))

    def test_effective_sample_size(self):
        """The effective number of samples should be sufficient."""
        self.run_check("effective_sample_size",
                       lambda posterior, prior, info: utils.check_effective_sample_size(len(posterior), info),
                       blocking=True)

    def test_multimodality(self):
        """Parameters should not have several modes (warning only)."""
        self.run_check("multimodality",
                       lambda posterior, prior, info: utils.check_multimodality(posterior))

    def test_derived_parameter_consistency(self):
        """Derived parameters should agree with the parameters they derive from."""
        self.run_check("parameter_consistency",
                       lambda posterior, prior, info: utils.check_derived_parameters(posterior),
                       blocking=True)

    def test_likelihood_sanity(self):
        """The log-likelihood and log-prior samples should be well behaved."""
        self.run_check("likelihood_sanity",
                       lambda posterior, prior, info: utils.check_likelihood(posterior),
                       blocking=True)

    def test_snr_plausibility(self):
        """The signal strength should be plausible for a near-threshold candidate (warning only)."""
        self.run_check("snr_plausibility",
                       lambda posterior, prior, info: utils.check_snr_plausibility(posterior, info))

    def test_sky_localization_quality(self):
        """Record the 90% sky area (informational)."""
        self.run_check("sky_localization",
                       lambda posterior, prior, info: utils.check_sky_localization(posterior))

    def test_network_snr_consistency(self):
        """The network SNR should be the quadrature sum of the detector SNRs."""
        self.run_check("network_snr",
                       lambda posterior, prior, info: utils.check_network_snr(posterior),
                       blocking=True)

    @classmethod
    def generate_event_html(cls, event_name, event_results, timestamp):
        """
        Generate HTML report for a single event showing only that event's results.

        Parameters
        ----------
        event_name : str
            Name of the event
        event_results : dict
            Dictionary of production results for this event
        timestamp : str
            Timestamp string for the report header

        Returns
        -------
        html : str
            Complete HTML document as a string
        """
        # Calculate summary statistics for this event only
        total_tests = 0
        status_counts = defaultdict(int)

        for production_results in event_results.values():
            for result in production_results.values():
                total_tests += 1
                status_counts[result['status']] += 1

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Automated PE Review Report - {event_name}</title>
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

        .header .event-title {{
            font-size: 1.5em;
            opacity: 0.95;
            margin-top: 10px;
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

        .production-section {{
            margin-bottom: 30px;
            border-left: 4px solid #667eea;
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
        <div class="event-title">{event_name}</div>
        <div class="timestamp">Generated: {timestamp}</div>
    </div>

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

    <div class="legend">
        <div class="legend-title">Test Descriptions</div>

        <div class="legend-item">
            <div class="legend-item-name">Spin Parameter Comparison</div>
            <div class="legend-item-desc">
                Compares spin parameter posteriors (a_1, a_2, chi_eff, chi_p) to their priors.
                Warns if posterior shows significantly more high-spin support (>0.7) than the
                prior, indicating data preference for high spins that may warrant attention.
            </div>
        </div>

        <div class="legend-item">
            <div class="legend-item-name">Sample Count Check</div>
            <div class="legend-item-desc">
                Verifies sufficient samples were collected. For nested sampling, samples are
                independent (IID). Thresholds: FAIL if <2000, WARNING if 2000-5000, PASS if ≥5000.
            </div>
        </div>

        <div class="legend-item">
            <div class="legend-item-name">Multimodality Check</div>
            <div class="legend-item-desc">
                Screens for multimodal posteriors by detecting large gaps in the distribution.
                Multimodality can indicate multiple solutions, poor convergence, or interesting
                astrophysics. Flagged cases warrant visual inspection of corner plots.
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

    <div class="event-section">
"""

        # Add results for each production in this event
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
</body>
</html>
"""

        return html

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
        # Generate a report per event that only shows that event's results
        for event_name, event_results in cls.test_results.items():
            # Create event-specific HTML report
            event_html = cls.generate_event_html(event_name, event_results, timestamp)

            for production_name in event_results.keys():
                report_path = get_report_path(event_name, production_name)

                # Create directory if it doesn't exist
                os.makedirs(os.path.dirname(report_path), exist_ok=True)

                # Write report
                try:
                    with open(report_path, 'w') as f:
                        f.write(event_html)
                    print(f"Generated automated review report: {report_path}")
                except Exception as e:
                    print(f"Warning: Could not write report to {report_path}: {e}")


if __name__ == '__main__':
    import unittest
    unittest.main()
