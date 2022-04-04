"""Module for testing Genomic Silent Mutation Validator."""
import unittest

from variation.validators import GenomicSilentMutation
from variation.classifiers import GenomicSilentMutationClassifier
from .validator_base import ValidatorBase


class TestGenomicSilentMutationValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Genomic Silent Mutation Validator."""

    def validator_instance(self):
        """Return genomic silent mutation instance."""
        return GenomicSilentMutation(*self.params)

    def classifier_instance(self):
        """Return the genomic silent mutation classifier instance."""
        return GenomicSilentMutationClassifier()

    def fixture_name(self):
        """Return the fixture name for genomic silent mutation."""
        return "genomic_silent_mutation"
