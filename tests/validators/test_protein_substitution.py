"""Module for testing Protein Substitution Validator."""
import unittest
from variant.validators import ProteinSubstitution
from variant.classifiers import ProteinSubstitutionClassifier
from .validator_base import ValidatorBase
from variant.data_sources import TranscriptMappings, SeqRepoAccess
from variant import PROJECT_ROOT


class TestProteinSubstitutionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Protein Substitution Validator."""

    def validator_instance(self):
        """Return protein substitution instance."""
        return ProteinSubstitution(SeqRepoAccess(
            f"{PROJECT_ROOT}/variant/data/seqrepo/latest"),
            TranscriptMappings(f"{PROJECT_ROOT}/variant/data"
                               f"/transcript_mapping.tsv")
        )

    def classifier_instance(self):
        """Return the protein substitution classifier instance."""
        return ProteinSubstitutionClassifier()

    def fixture_name(self):
        """Return the fixture name for protein substitution."""
        return 'protein_substitution'
