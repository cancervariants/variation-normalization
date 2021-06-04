"""Module for testing Coding DNA Substitution Translator."""
import unittest
from variation.classifiers import CodingDNASubstitutionClassifier
from variation.translators import CodingDNASubstitution
from variation.validators import CodingDNASubstitution as CDNASUB_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache
from variation.data_sources import SeqRepoAccess, TranscriptMappings
from variation import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestCodingDNASubstitutionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Coding DNA Substitution Translator."""

    def classifier_instance(self):
        """Return coding DNA substitution instance."""
        return CodingDNASubstitutionClassifier()

    def validator_instance(self):
        """Return coding DNA substitution instance."""
        return CDNASUB_V(SeqRepoAccess(SEQREPO_DATA_PATH),
                         TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),
                         GeneSymbol(GeneSymbolCache())
                         )

    def translator_instance(self):
        """Return coding DNA substitution instance."""
        return CodingDNASubstitution()

    def fixture_name(self):
        """Return the fixture name for coding DNA substitution."""
        return 'coding_dna_substitution'
