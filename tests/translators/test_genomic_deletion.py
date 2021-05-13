"""Module for testing Genomic Deletion Translator."""
import unittest
from variant.classifiers import GenomicDeletionClassifier
from variant.translators import GenomicDeletion
from variant.validators import GenomicDeletion as GD_V
from .translator_base import TranslatorBase
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache
from variant.data_sources import SeqRepoAccess, TranscriptMappings
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestGenomicDeletionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Genomic Deletion Translator."""

    def classifier_instance(self):
        """Return genomic deletion instance."""
        return GenomicDeletionClassifier()

    def validator_instance(self):
        """Return coding DNA delins instance."""
        return GD_V(SeqRepoAccess(SEQREPO_DATA_PATH),
                    TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),
                    GeneSymbol(GeneSymbolCache())
                    )

    def translator_instance(self):
        """Return genomic deletion instance."""
        return GenomicDeletion()

    def fixture_name(self):
        """Return the fixture name for genomic deletion."""
        return 'genomic_deletion'
