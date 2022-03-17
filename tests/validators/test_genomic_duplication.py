"""Module for testing Genomic Duplication Validator."""
import unittest

from variation.validators import GenomicDuplication
from variation.classifiers import GenomicDuplicationClassifier
from .validator_base import ValidatorBase


class TestGenomicDuplicationValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Genomic Duplication Validator."""

    def validator_instance(self):
        """Return genomic duplication instance."""
        return GenomicDuplication(*self.params)

    def classifier_instance(self):
        """Return the genomic duplication classifier instance."""
        return GenomicDuplicationClassifier()

    def fixture_name(self):
        """Return the fixture name for genomic duplication."""
        return "genomic_duplication"
