"""Module for testing Polypeptide Truncation Validator."""
import unittest
from variant.validators import PolypeptideTruncation
from variant.classifiers import PolypeptideTruncationClassifier
from .validator_base import ValidatorBase
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from variant.data_sources import TranscriptMappings, SeqRepoAccess
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestPolypeptideTruncationValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Polypeptide Truncation Validator."""

    def validator_instance(self):
        """Return Polypeptide Truncation instance."""
        return PolypeptideTruncation(SeqRepoAccess(SEQREPO_DATA_PATH),
                                     TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),  # noqa: E501
                                     GeneSymbol(GeneSymbolCache()),
                                     AminoAcidCache())

    def classifier_instance(self):
        """Return the Polypeptide Truncation classifier instance."""
        return PolypeptideTruncationClassifier()

    def fixture_name(self):
        """Return the fixture name for Polypeptide Truncation."""
        return 'polypeptide_truncation'
