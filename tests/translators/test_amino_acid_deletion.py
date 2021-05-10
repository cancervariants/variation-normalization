"""Module for testing amino acid deletion Translator."""
import unittest
from variant.classifiers import AminoAcidDeletionClassifier
from variant.translators import AminoAcidDeletion
from variant.validators import AminoAcidDeletion as AAD_V
from .translator_base import TranslatorBase
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from variant.data_sources import SeqRepoAccess, TranscriptMappings
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestAminoAcidDeletionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the amino acid deletion Translator."""

    def classifier_instance(self):
        """Return amino acid deletion instance."""
        return AminoAcidDeletionClassifier()

    def validator_instance(self):
        """Return amino acid deletion instance."""
        return AAD_V(SeqRepoAccess(SEQREPO_DATA_PATH),
                     TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),
                     GeneSymbol(GeneSymbolCache()),
                     AminoAcidCache()
                     )

    def translator_instance(self):
        """Return amino acid deletion instance."""
        return AminoAcidDeletion()

    def fixture_name(self):
        """Return the fixture name for amino acid deletion."""
        return 'amino_acid_deletion'
