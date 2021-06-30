"""Module for testing Silent Mutation Validator."""
import unittest
from variation.validators import SilentMutation
from variation.classifiers import SilentMutationClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript


class TestSilentMutationValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Silent Mutation Validator."""

    def validator_instance(self):
        """Return Silent Mutation instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        return SilentMutation(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), UTA()),
            AminoAcidCache())

    def classifier_instance(self):
        """Return the Silent Mutation classifier instance."""
        return SilentMutationClassifier()

    def fixture_name(self):
        """Return the fixture name for Silent Mutation."""
        return 'silent_mutation'
