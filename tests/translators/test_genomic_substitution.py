"""Module for testing Genomic Substitution Translator."""
import unittest
from variation.classifiers import GenomicSubstitutionClassifier
from variation.translators import GenomicSubstitution
from variation.validators import GenomicSubstitution as GSUB_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript


class TestGenomicSubstitutionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Genomic Substitution Translator."""

    def classifier_instance(self):
        """Return genomic substitution instance."""
        return GenomicSubstitutionClassifier()

    def validator_instance(self):
        """Return genomic substitution instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        return GSUB_V(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta
        )

    def translator_instance(self):
        """Return genomic substitution instance."""
        return GenomicSubstitution()

    def fixture_name(self):
        """Return the fixture name for genomic substitution."""
        return 'genomic_substitution'
