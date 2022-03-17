"""Module for testing Genomic Insertion Validator."""
import unittest

from variation.validators import GenomicInsertion
from variation.classifiers import GenomicInsertionClassifier
from .validator_base import ValidatorBase


class TestGenomicInsertionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the GenomicInsertion Validator."""

    def validator_instance(self):
        """Return genomic insertion instance."""
        return GenomicInsertion(*self.params)

    def classifier_instance(self):
        """Return the genomic insertion classifier instance."""
        return GenomicInsertionClassifier()

    def fixture_name(self):
        """Return the fixture name for genomic insertion."""
        return "genomic_insertion"
