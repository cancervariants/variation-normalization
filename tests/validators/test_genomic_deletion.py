"""Module for testing Genomic Deletion Validator."""
import unittest

from variation.validators import GenomicDeletion
from variation.classifiers import GenomicDeletionClassifier
from .validator_base import ValidatorBase


class TestGenomicDeletionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the GenomicDeletion Validator."""

    def validator_instance(self):
        """Return genomic deletion instance."""
        return GenomicDeletion(*self.params)

    def classifier_instance(self):
        """Return the genomic deletion classifier instance."""
        return GenomicDeletionClassifier()

    def fixture_name(self):
        """Return the fixture name for genomic deletion."""
        return "genomic_deletion"
