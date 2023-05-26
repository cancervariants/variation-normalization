"""Module for testing Genomic Reference Agree Validator."""
import unittest

from variation.validators import GenomicReferenceAgree
from variation.classifiers import GenomicReferenceAgreeClassifier
from .validator_base import ValidatorBase


class TestGenomicReferenceAgreeValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Genomic Reference Agree Validator."""

    def validator_instance(self):
        """Return genomic reference agree instance."""
        return GenomicReferenceAgree(*self.params)

    def classifier_instance(self):
        """Return the genomic reference agree classifier instance."""
        return GenomicReferenceAgreeClassifier()

    def fixture_name(self):
        """Return the fixture name for genomic reference agree."""
        return "genomic_reference_agree"
