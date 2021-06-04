"""Module for testing Genomic Substitution Translator."""
import unittest
from variation.classifiers import GenomicSubstitutionClassifier
from variation.translators import GenomicSubstitution
from variation.validators import GenomicSubstitution as GSUB_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache
from variation.data_sources import SeqRepoAccess, TranscriptMappings
from variation import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestGenomicSubstitutionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Genomic Substitution Translator."""

    def classifier_instance(self):
        """Return genomic substitution instance."""
        return GenomicSubstitutionClassifier()

    def validator_instance(self):
        """Return genomic substitution instance."""
        return GSUB_V(SeqRepoAccess(SEQREPO_DATA_PATH),
                      TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),
                      GeneSymbol(GeneSymbolCache())
                      )

    def translator_instance(self):
        """Return genomic substitution instance."""
        return GenomicSubstitution()

    def fixture_name(self):
        """Return the fixture name for genomic substitution."""
        return 'genomic_substitution'
