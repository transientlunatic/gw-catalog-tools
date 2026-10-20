"""
The tests which run in an asimov project.

Importing this module needs asimov, and, as asimov does when it is imported,
reads the ledger of the project in the current directory.

To run the tests, put them in the project's tests directory with
``catalog-tools install-tests``, or import them into a test module of your own::

    from catalog_tools.asimov_tests import *
"""

from .automated_review import TestAutomatedReview
from .samples import TestSampleLoss
from .waveform import TestWaveformConfiguration

__all__ = ["TestAutomatedReview", "TestSampleLoss", "TestWaveformConfiguration"]
