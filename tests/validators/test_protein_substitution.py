"""Module for testing Protein Substitution Validator."""
import unittest
from variant.validators import ProteinSubstitution
from variant.classifiers import ProteinSubstitutionClassifier
from .validator_base import ValidatorBase
from variant.data_sources import TranscriptMappings, SeqRepoAccess
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestProteinSubstitutionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Protein Substitution Validator."""

    def validator_instance(self):
        """Return protein substitution instance."""
        return ProteinSubstitution(SeqRepoAccess(SEQREPO_DATA_PATH),
                                   TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH)
                                   )

    def classifier_instance(self):
        """Return the protein substitution classifier instance."""
        return ProteinSubstitutionClassifier()

    def fixture_name(self):
        """Return the fixture name for protein substitution."""
        return 'protein_substitution'
