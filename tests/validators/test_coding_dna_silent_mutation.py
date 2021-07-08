"""Module for testing Coding DNA Silent Mutation Validator."""
import unittest
from variation.validators import CodingDNASilentMutation
from variation.classifiers import CodingDNASilentMutationClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript


class TestCodingDNASilentMutationValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Coding DNA Silent Mutation Validator."""

    def validator_instance(self):
        """Return coding DNA silent mutation instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        return CodingDNASilentMutation(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta
        )

    def classifier_instance(self):
        """Return the coding DNA silent mutation classifier instance."""
        return CodingDNASilentMutationClassifier()

    def fixture_name(self):
        """Return the fixture name for coding DNA silent mutation."""
        return 'coding_dna_silent_mutation'
