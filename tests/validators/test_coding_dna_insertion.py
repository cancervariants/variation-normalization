"""Module for testing Coding DNA Insertion Validator."""
import unittest
from variation.validators import CodingDNAInsertion
from variation.classifiers import CodingDNAInsertionClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript


class TestCodingDNAInsertionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the CodingDNAInsertion Validator."""

    def validator_instance(self):
        """Return coding dna insertion instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        return CodingDNAInsertion(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), UTA())
        )

    def classifier_instance(self):
        """Return the coding dna insertion classifier instance."""
        return CodingDNAInsertionClassifier()

    def fixture_name(self):
        """Return the fixture name for coding dna insertion."""
        return 'coding_dna_insertion'
