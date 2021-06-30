"""Module for testing amino acid insertion Translator."""
import unittest
from variation.classifiers import AminoAcidInsertionClassifier
from variation.translators import AminoAcidInsertion
from variation.validators import AminoAcidInsertion as AAI_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript


class TestAminoAcidInsertionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the amino acid insertion Translator."""

    def classifier_instance(self):
        """Return amino acid insertion instance."""
        return AminoAcidInsertionClassifier()

    def validator_instance(self):
        """Return amino acid insertion instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        return AAI_V(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), UTA()),
            AminoAcidCache())

    def translator_instance(self):
        """Return amino acid insertion instance."""
        return AminoAcidInsertion()

    def fixture_name(self):
        """Return the fixture name for amino acid insertion."""
        return 'amino_acid_insertion'
