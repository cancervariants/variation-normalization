"""Module for testing Genomic Uncertain Deletion Validator."""
import unittest

from variation.validators import GenomicUncertainDeletion
from variation.classifiers import GenomicUncertainDeletionClassifier
from .validator_base import ValidatorBase


class TestGenomicUncertainDeletionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Genomic Uncertain Deletion Validator."""

    def validator_instance(self):
        """Return genomic uncertain deletion instance."""
        return GenomicUncertainDeletion(*self.params)

    def classifier_instance(self):
        """Return the genomic uncertain deletion classifier instance."""
        return GenomicUncertainDeletionClassifier()

    def fixture_name(self):
        """Return the fixture name for genomic uncertain deletion."""
        return "genomic_uncertain_deletion"
