"""Module for testing Silent Mutation Validator."""
import unittest
from variant.validators import SilentMutation
from variant.classifiers import SilentMutationClassifier
from .validator_base import ValidatorBase
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from variant.data_sources import TranscriptMappings, SeqRepoAccess
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestSilentMutationValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Silent Mutation Validator."""

    def validator_instance(self):
        """Return Silent Mutation instance."""
        return SilentMutation(SeqRepoAccess(SEQREPO_DATA_PATH),
                              TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),  # noqa: E501
                              GeneSymbol(GeneSymbolCache()),
                              AminoAcidCache())

    def classifier_instance(self):
        """Return the Silent Mutation classifier instance."""
        return SilentMutationClassifier()

    def fixture_name(self):
        """Return the fixture name for Silent Mutation."""
        return 'silent_mutation'
