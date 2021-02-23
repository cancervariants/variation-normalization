"""Module for testing Amino Acid Substitution Validator."""
import unittest
from variant.validators import PolypeptideTruncation
from variant.classifiers import PolypeptideTruncationClassifier
from .validator_base import ValidatorBase
from variant.data_sources import TranscriptMappings, SeqRepoAccess
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestPolypeptideTruncationValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Protein Substitution Validator."""

    def validator_instance(self):
        """Return amino acid substitution instance."""
        return PolypeptideTruncation(SeqRepoAccess(SEQREPO_DATA_PATH),
                                     TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH))  # noqa: E501

    def classifier_instance(self):
        """Return the protein substitution classifier instance."""
        return PolypeptideTruncationClassifier()

    def fixture_name(self):
        """Return the fixture name for protein substitution."""
        return 'polypeptide_truncation'
