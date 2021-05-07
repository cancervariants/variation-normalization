"""Module for testing amino acid delins Translator."""
import unittest
from variant.classifiers import AminoAcidDelInsClassifier
from variant.translators import AminoAcidDelIns
from variant.validators import AminoAcidDelIns as AAD_V
from .translator_base import TranslatorBase
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from variant.data_sources import SeqRepoAccess, TranscriptMappings
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestAminoAcidDelInsTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the amino acid delins Translator."""

    def classifier_instance(self):
        """Return amino acid delins instance."""
        return AminoAcidDelInsClassifier()

    def validator_instance(self):
        """Return amino acid delins instance."""
        return AAD_V(SeqRepoAccess(SEQREPO_DATA_PATH),
                     TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),
                     GeneSymbol(GeneSymbolCache()),
                     AminoAcidCache()
                     )

    def translator_instance(self):
        """Return amino acid delins instance."""
        return AminoAcidDelIns()

    def fixture_name(self):
        """Return the fixture name for amino acid delins."""
        return 'amino_acid_delins'
