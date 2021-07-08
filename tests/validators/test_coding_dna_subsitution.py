"""Module for testing Coding DNA Substitution Validator."""
import unittest
from variation.validators import CodingDNASubstitution
from variation.classifiers import CodingDNASubstitutionClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript


class TestCodingDNASubstitutionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Coding DNA Substitution Validator."""

    def validator_instance(self):
        """Return coding DNA substitution instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        return CodingDNASubstitution(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta
        )

    def classifier_instance(self):
        """Return the coding DNA substitution classifier instance."""
        return CodingDNASubstitutionClassifier()

    def fixture_name(self):
        """Return the fixture name for coding DNA substitution."""
        return 'coding_dna_substitution'
