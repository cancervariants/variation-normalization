"""Module for testing amino acid insertion Translator."""
import unittest
from variant.classifiers import AminoAcidInsertionClassifier
from variant.translators import AminoAcidInsertion
from variant.validators import AminoAcidInsertion as AAI_V
from .translator_base import TranslatorBase
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from variant.data_sources import SeqRepoAccess, TranscriptMappings
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestAminoAcidInsertionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the amino acid insertion Translator."""

    def classifier_instance(self):
        """Return amino acid insertion instance."""
        return AminoAcidInsertionClassifier()

    def validator_instance(self):
        """Return amino acid insertion instance."""
        return AAI_V(SeqRepoAccess(SEQREPO_DATA_PATH),
                     TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),
                     GeneSymbol(GeneSymbolCache()),
                     AminoAcidCache()
                     )

    def translator_instance(self):
        """Return amino acid insertion instance."""
        return AminoAcidInsertion()

    def fixture_name(self):
        """Return the fixture name for amino acid insertion."""
        return 'amino_acid_insertion'
