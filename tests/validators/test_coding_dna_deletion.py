"""Module for testing Coding DNA Deletion Validator."""
import unittest
from variation.validators import CodingDNADeletion
from variation.classifiers import CodingDNADeletionClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript


class TestCodingDNADeletionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the CodingDNADeletion Validator."""

    def validator_instance(self):
        """Return coding dna deletion instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        return CodingDNADeletion(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta
        )

    def classifier_instance(self):
        """Return the coding dna deletion classifier instance."""
        return CodingDNADeletionClassifier()

    def fixture_name(self):
        """Return the fixture name for coding dna deletion."""
        return 'coding_dna_deletion'
