"""Module for testing Genomic DelIns Validator."""
import unittest

from variation.validators import GenomicDelIns
from variation.classifiers import GenomicDelInsClassifier
from .validator_base import ValidatorBase


class TestGenomicDelInsValidator(ValidatorBase, unittest.TestCase):
    """A class to test the GenomicDelIns Validator."""

    def validator_instance(self):
        """Return genomic delins instance."""
        return GenomicDelIns(*self.params)

    def classifier_instance(self):
        """Return the genomic delins classifier instance."""
        return GenomicDelInsClassifier()

    def fixture_name(self):
        """Return the fixture name for genomic delins."""
        return "genomic_delins"
