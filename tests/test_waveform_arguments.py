from asimov.testing import AsimovTest
import configparser
from bilby_pipe import bilbyargparser
import ast
#import bilbyargparser.BilbyConfigFileParser as bilby_parse
import os
import glob


class TestBilbyAnalyses(AsimovTest):
    def test_spin_taylor_args(self):
        """Verify the component chirp mass prior is correct."""

        for event in self.events:
            for production in event.productions:
                if (production.pipeline.name == "Bilby") \
                   and (production.status_str.lower() == "uploaded") \
                   and ("IMRPhenomXPHM-SpinTaylor" in production.name)\
                   and (production.review.status == "APPROVED"):
                    with self.subTest(event=event.name, production=production.name):
                        # Read the config file
                        parser = bilbyargparser.BilbyConfigFileParser()
                        with open(f'{production.meta["working directory"]}/{production.name}/{production.name}_config_complete.ini', "r") as filename:
                            config = parser.parse(filename)[0]
                            self.assertTrue("'PhenomXHMReleaseVersion': 122022" in config['waveform-arguments-dict'])
                            self.assertTrue("'PhenomXPFinalSpinMod': 2" in config['waveform-arguments-dict'])
                            self.assertTrue("'PhenomXPrecVersion': 320" in config['waveform-arguments-dict'])

    def test_l_max(self):
        """Test for issues with l-max."""

        for event in self.events:
            for production in event.productions:
                if (production.pipeline.name == "Bilby") \
                   and (production.status_str.lower() == "uploaded") \
                   and ("SEOBNRv5PHM" in production.name)\
                   and (production.review.status == "APPROVED"):
                    with self.subTest(event=event.name, production=production.name):

                        if (production.meta['waveform']['arguments']['lmax_nyquist'] != 1):
                            err_files = glob.glob(f'{production.meta["working directory"]}/{production.name}/log_data_analysis/*.err')
                            for error_file in err_files:
                                with open(error_file, "r") as error_file_handle:
                                    data = error_file_handle.read()
                                    count = data.count("Ringdown frequency of (3,3) mode greater than maximum frequency from Nyquist theorem")
                                    if count > 0:
                                        print(event, count)
                                    self.assertTrue(count < 5)                            

    # def test_xo4q_env(self):
    #     """Verify the environment used for the XO4a waveform is correct."""

    #     for event in self.events:
    #         for production in event.productions:
    #             if (production.pipeline.name.lower() == "bilby") \
    #                and (production.status_str.lower() == "uploaded") \
    #                and ("IMRPhenomXO4a" in production.name)\
    #                and (production.review.status == "APPROVED"):
    #                 with self.subTest(event=event.name, production=production.name):
    #                     # Read the config file
    #                     parser = bilbyargparser.BilbyConfigFileParser()
    #                     with open(f'{production.meta["working directory"]}/{production.name}/{production.name}_config_complete.ini', "r") as filename:
    #                         config = parser.parse(filename)[0]
    #                         self.assertTrue("/cvmfs/software.igwn.org/conda/envs/igwn-py310-20250129" in config['analysis-executable'])
                        
if __name__ == '__main__':
    unittest.main()
