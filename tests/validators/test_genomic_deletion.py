"""Module for testing Genomic Deletion Validator."""
import unittest
from variant.validators import GenomicDeletion
from variant.classifiers import GenomicDeletionClassifier
from .validator_base import ValidatorBase
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache
from variant.data_sources import TranscriptMappings, SeqRepoAccess
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestGenomicDeletionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the GenomicDeletion Validator."""

    def validator_instance(self):
        """Return genomic deletion instance."""
        return GenomicDeletion(SeqRepoAccess(SEQREPO_DATA_PATH),
                               TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),  # noqa: E501
                               GeneSymbol(GeneSymbolCache()))

    def classifier_instance(self):
        """Return the genomic deletion classifier instance."""
        return GenomicDeletionClassifier()

    def fixture_name(self):
        """Return the fixture name for genomic deletion."""
        return 'genomic_deletion'
