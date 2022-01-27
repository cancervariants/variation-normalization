"""Module for testing Amino Acid Substitution Validator."""
import unittest
from variation.validators import AminoAcidSubstitution
from variation.classifiers import AminoAcidSubstitutionClassifier
from .validator_base import ValidatorBase


class TestAminoAcidSubstitutionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Protein Substitution Validator."""

    def validator_instance(self):
        """Return amino acid substitution instance."""
        return AminoAcidSubstitution(*self.aa_params)

    def classifier_instance(self):
        """Return the protein substitution classifier instance."""
        return AminoAcidSubstitutionClassifier()

    def fixture_name(self):
        """Return the fixture name for protein substitution."""
        return 'amino_acid_substitution'
