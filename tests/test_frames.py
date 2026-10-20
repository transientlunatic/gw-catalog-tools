import unittest

from asimov.testing import AsimovTest
from bilby_pipe import bilbyargparser

from review_utils import is_bilby, is_uploaded

class TestBilbyAnalyses(AsimovTest):
    def test_phenomxhm_release_version(self):
        """Verify the IMRPhenomXPHM release version in each analysis' configuration."""

        for event in self.events:
            for production in event.productions:
                if is_bilby(production) and is_uploaded(production) and ("IMRPhenomXPHM-SpinTaylor" in production.name):
                    with self.subTest(event=event.name, production=production.name):
                        # Read the config file
                        parser = bilbyargparser.BilbyConfigFileParser()
                        with open(f'{production.meta["working directory"]}/{production.name}/{production.name}_config_complete.ini', "r") as filename:
                            config = parser.parse(filename)[0]
                            self.assertTrue("'PhenomXHMReleaseVersion': 122022" in config['waveform-arguments-dict'],
                                            msg=f"{event.name} {production.name}")

                            
if __name__ == '__main__':
    unittest.main()
