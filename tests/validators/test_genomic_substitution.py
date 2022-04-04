"""Module for testing Genomic Substitution Validator."""
import unittest

from variation.validators import GenomicSubstitution
from variation.classifiers import GenomicSubstitutionClassifier
from .validator_base import ValidatorBase


class TestGenomicSubstitutionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Genomic Substitution Validator."""

    def validator_instance(self):
        """Return genomic substitution instance."""
        return GenomicSubstitution(*self.params)

    def classifier_instance(self):
        """Return the genomic substitution classifier instance."""
        return GenomicSubstitutionClassifier()

    def fixture_name(self):
        """Return the fixture name for genomic substitution."""
        return "genomic_substitution"
