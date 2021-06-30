"""Module for testing Coding DNA Substitution Translator."""
import unittest
from variation.classifiers import CodingDNASubstitutionClassifier
from variation.translators import CodingDNASubstitution
from variation.validators import CodingDNASubstitution as CDNASUB_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript


class TestCodingDNASubstitutionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Coding DNA Substitution Translator."""

    def classifier_instance(self):
        """Return coding DNA substitution instance."""
        return CodingDNASubstitutionClassifier()

    def validator_instance(self):
        """Return coding DNA substitution instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        return CDNASUB_V(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), UTA())
        )

    def translator_instance(self):
        """Return coding DNA substitution instance."""
        return CodingDNASubstitution()

    def fixture_name(self):
        """Return the fixture name for coding DNA substitution."""
        return 'coding_dna_substitution'
