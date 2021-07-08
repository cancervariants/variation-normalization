"""Module for testing Genomic Deletion Translator."""
import unittest
from variation.classifiers import GenomicDeletionClassifier
from variation.translators import GenomicDeletion
from variation.validators import GenomicDeletion as GD_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript


class TestGenomicDeletionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Genomic Deletion Translator."""

    def classifier_instance(self):
        """Return genomic deletion instance."""
        return GenomicDeletionClassifier()

    def validator_instance(self):
        """Return coding DNA deletion instance."""
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
        """Return genomic deletion instance."""
        return GenomicDeletion()

    def fixture_name(self):
        """Return the fixture name for genomic deletion."""
        return 'genomic_deletion'
