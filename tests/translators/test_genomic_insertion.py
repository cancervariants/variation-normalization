"""Module for testing Genomic Insertion Translator."""
import unittest
from variation.classifiers import GenomicInsertionClassifier
from variation.translators import GenomicInsertion
from variation.validators import GenomicInsertion as GD_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache
from variation.data_sources import SeqRepoAccess, TranscriptMappings
from variation import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestGenomicInsertionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Genomic Insertion Translator."""

    def classifier_instance(self):
        """Return genomic insertion instance."""
        return GenomicInsertionClassifier()

    def validator_instance(self):
        """Return coding DNA insertion instance."""
        return GD_V(SeqRepoAccess(SEQREPO_DATA_PATH),
                    TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),
                    GeneSymbol(GeneSymbolCache())
                    )

    def translator_instance(self):
        """Return genomic insertion instance."""
        return GenomicInsertion()

    def fixture_name(self):
        """Return the fixture name for genomic insertion."""
        return 'genomic_insertion'
