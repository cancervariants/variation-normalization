"""Module for testing Coding DNA Silent Mutation Validator."""
import unittest
from variant.validators import CodingDNASilentMutation
from variant.classifiers import CodingDNASilentMutationClassifier
from .validator_base import ValidatorBase
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache
from variant.data_sources import TranscriptMappings, SeqRepoAccess
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestCodingDNASilentMutationValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Coding DNA Silent Mutation Validator."""

    def validator_instance(self):
        """Return coding DNA silent mutation instance."""
        return CodingDNASilentMutation(SeqRepoAccess(SEQREPO_DATA_PATH),
                                       TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),  # noqa: E501
                                       GeneSymbol(GeneSymbolCache()))

    def classifier_instance(self):
        """Return the coding DNA silent mutation classifier instance."""
        return CodingDNASilentMutationClassifier()

    def fixture_name(self):
        """Return the fixture name for coding DNA silent mutation."""
        return 'coding_dna_silent_mutation'
