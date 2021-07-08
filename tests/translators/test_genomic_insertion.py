"""Module for testing Genomic Insertion Translator."""
import unittest
from variation.classifiers import GenomicInsertionClassifier
from variation.translators import GenomicInsertion
from variation.validators import GenomicInsertion as GD_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript


class TestGenomicInsertionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Genomic Insertion Translator."""

    def classifier_instance(self):
        """Return genomic insertion instance."""
        return GenomicInsertionClassifier()

    def validator_instance(self):
        """Return coding DNA insertion instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        return GD_V(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta
        )

    def translator_instance(self):
        """Return genomic insertion instance."""
        return GenomicInsertion()

    def fixture_name(self):
        """Return the fixture name for genomic insertion."""
        return 'genomic_insertion'
