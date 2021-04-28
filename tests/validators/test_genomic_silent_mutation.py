"""Module for testing Genomic Silent Mutation Validator."""
import unittest
from variant.validators import GenomicSilentMutation
from variant.classifiers import GenomicSilentMutationClassifier
from .validator_base import ValidatorBase
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache
from variant.data_sources import TranscriptMappings, SeqRepoAccess
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestGenomicSilentMutationValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Genomic Silent Mutation Validator."""

    def validator_instance(self):
        """Return genomic silent mutation instance."""
        return GenomicSilentMutation(SeqRepoAccess(SEQREPO_DATA_PATH),
                                     TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),  # noqa: E501
                                     GeneSymbol(GeneSymbolCache()))

    def classifier_instance(self):
        """Return the genomic silent mutation classifier instance."""
        return GenomicSilentMutationClassifier()

    def fixture_name(self):
        """Return the fixture name for genomic silent mutation."""
        return 'genomic_silent_mutation'
