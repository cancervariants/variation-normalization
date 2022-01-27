"""Module for testing Amino Acid Insertion Validator."""
import unittest
from variation.validators import AminoAcidInsertion
from variation.classifiers import AminoAcidInsertionClassifier
from .validator_base import ValidatorBase


class TestAminoAcidInsertionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Amino Acid Insertion Validator."""

    def validator_instance(self):
        """Return amino acid insertion instance."""
        return AminoAcidInsertion(*self.aa_params)

    def classifier_instance(self):
        """Return the amino acid insertion classifier instance."""
        return AminoAcidInsertionClassifier()

    def fixture_name(self):
        """Return the fixture name for amino acid insertion."""
        return 'amino_acid_insertion'
