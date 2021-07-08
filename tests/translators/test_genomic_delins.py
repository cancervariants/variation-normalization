"""Module for testing Genomic DelIns Translator."""
import unittest
from variation.classifiers import GenomicDelInsClassifier
from variation.translators import GenomicDelIns
from variation.validators import GenomicDelIns as GENOMICDELINS_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript


class TestGenomicDelInsTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Genomic DelIns Translator."""

    def classifier_instance(self):
        """Return Genomic delins instance."""
        return GenomicDelInsClassifier()

    def validator_instance(self):
        """Return genomic delins instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        return GENOMICDELINS_V(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta
        )

    def translator_instance(self):
        """Return genomic delins instance."""
        return GenomicDelIns()

    def fixture_name(self):
        """Return the fixture name for genomic delins."""
        return 'genomic_delins'
