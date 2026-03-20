from asimov.testing import AsimovTest
from bilby_pipe import bilbyargparser

class TestBilbyAnalyses(AsimovTest):
    def test_component_chirp_mass(self):
        """Verify the component chirp mass prior is correct."""

        for event in self.events:
            for production in event.productions:
                if (production.pipeline.name == "Bilby") and (production.status_str.lower() == "uploaded") and ("IMRPhenomXPHM-SpinTaylor" in production.name):
                    with self.subTest(event=event.name, production=production.name):
                        # Read the config file
                        parser = bilbyargparser.BilbyConfigFileParser()
                        with open(f'{production.meta["working directory"]}/{production.name}/{production.name}_config_complete.ini', "r") as filename:
                            config = parser.parse(filename)[0]
                            self.assertTrue("'PhenomXHMReleaseVersion': 122022" in config['waveform-arguments-dict'])


