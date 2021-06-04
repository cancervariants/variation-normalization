"""Module for testing polypeptide truncation Translator."""
import unittest
from variation.classifiers import PolypeptideTruncationClassifier
from variation.translators import PolypeptideTruncation
from variation.validators import PolypeptideTruncation as PT_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from variation.data_sources import SeqRepoAccess, TranscriptMappings
from variation import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestPolypeptideTruncationTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the polypeptide truncation Translator."""

    def classifier_instance(self):
        """Return polypeptide truncation instance."""
        return PolypeptideTruncationClassifier()

    def validator_instance(self):
        """Return polypeptide truncation instance."""
        return PT_V(SeqRepoAccess(SEQREPO_DATA_PATH),
                    TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),
                    GeneSymbol(GeneSymbolCache()),
                    AminoAcidCache()
                    )

    def translator_instance(self):
        """Return polypeptide truncation instance."""
        return PolypeptideTruncation()

    def fixture_name(self):
        """Return the fixture name for polypeptide truncation."""
        return 'polypeptide_truncation'
