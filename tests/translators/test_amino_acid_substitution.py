"""Module for testing amino acid Substitution Translator."""
import unittest
from variation.classifiers import AminoAcidSubstitutionClassifier
from variation.translators import AminoAcidSubstitution
from variation.validators import AminoAcidSubstitution as AASUB_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript


class TestAminoAcidSubstitutionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the amino acid Substitution Translator."""

    def classifier_instance(self):
        """Return amino acid substitution instance."""
        return AminoAcidSubstitutionClassifier()

    def validator_instance(self):
        """Return amino acid substitution instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        return AASUB_V(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), UTA()),
            AminoAcidCache())

    def translator_instance(self):
        """Return amino acid substitution instance."""
        return AminoAcidSubstitution()

    def fixture_name(self):
        """Return the fixture name for amino acid substitution."""
        return 'amino_acid_substitution'
