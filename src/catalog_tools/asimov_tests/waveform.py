"""Check the waveform configuration of the analyses."""

import glob

from asimov.testing import AsimovTest
from bilby_pipe import bilbyargparser

from catalog_tools.review import is_approved, is_bilby, is_uploaded


def read_config(production):
    """The complete bilby configuration of an analysis, as a dictionary."""
    parser = bilbyargparser.BilbyConfigFileParser()
    path = f'{production.meta["working directory"]}/{production.name}/{production.name}_config_complete.ini'
    with open(path, "r") as handle:
        return parser.parse(handle)[0]


class TestWaveformConfiguration(AsimovTest):
    def test_phenomxhm_release_version(self):
        """Verify the IMRPhenomXPHM release version in each analysis' configuration."""
        for event in self.events:
            for production in event.productions:
                if is_bilby(production) and is_uploaded(production) \
                   and "IMRPhenomXPHM-SpinTaylor" in production.name:
                    with self.subTest(event=event.name, production=production.name):
                        config = read_config(production)
                        self.assertTrue("'PhenomXHMReleaseVersion': 122022" in config['waveform-arguments-dict'],
                                        msg=f"{event.name} {production.name}")

    def test_spin_taylor_args(self):
        """Verify the IMRPhenomXPHM-SpinTaylor waveform arguments of approved analyses."""
        for event in self.events:
            for production in event.productions:
                if is_bilby(production) and is_uploaded(production) \
                   and "IMRPhenomXPHM-SpinTaylor" in production.name \
                   and is_approved(production):
                    with self.subTest(event=event.name, production=production.name):
                        config = read_config(production)
                        for argument in ("'PhenomXHMReleaseVersion': 122022",
                                         "'PhenomXPFinalSpinMod': 2",
                                         "'PhenomXPrecVersion': 320"):
                            self.assertTrue(argument in config['waveform-arguments-dict'],
                                            msg=f"{event.name} {production.name}: missing {argument}")

    def test_l_max(self):
        """Test for issues with l-max in approved SEOBNRv5PHM analyses."""
        for event in self.events:
            for production in event.productions:
                if is_bilby(production) and is_uploaded(production) \
                   and "SEOBNRv5PHM" in production.name \
                   and is_approved(production):
                    with self.subTest(event=event.name, production=production.name):
                        if production.meta['waveform']['arguments']['lmax_nyquist'] != 1:
                            error_files = glob.glob(
                                f'{production.meta["working directory"]}/{production.name}/log_data_analysis/*.err')
                            for error_file in error_files:
                                with open(error_file, "r") as handle:
                                    count = handle.read().count(
                                        "Ringdown frequency of (3,3) mode greater than maximum "
                                        "frequency from Nyquist theorem")
                                self.assertTrue(count < 5, msg=f"{event.name} {production.name}")
