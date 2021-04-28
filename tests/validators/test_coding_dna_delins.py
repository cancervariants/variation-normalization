"""Module for testing Coding DNA DelIns Validator."""
import unittest
from variant.validators import CodingDNADelIns
from variant.classifiers import CodingDNADelInsClassifier
from .validator_base import ValidatorBase
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache
from variant.data_sources import TranscriptMappings, SeqRepoAccess
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestCodingDNADelInsValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Coding DNA DelIns Validator."""

    def validator_instance(self):
        """Return coding DNA delins instance."""
        return CodingDNADelIns(SeqRepoAccess(SEQREPO_DATA_PATH),
                               TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),  # noqa: E501
                               GeneSymbol(GeneSymbolCache()))

    def classifier_instance(self):
        """Return the coding DNA delins classifier instance."""
        return CodingDNADelInsClassifier()

    def fixture_name(self):
        """Return the fixture name for coding DNA delins."""
        return 'coding_dna_delins'
