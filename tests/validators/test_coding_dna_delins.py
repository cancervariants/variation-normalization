"""Module for testing Coding DNA DelIns Validator."""
import unittest
from variation.validators import CodingDNADelIns
from variation.classifiers import CodingDNADelInsClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript


class TestCodingDNADelInsValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Coding DNA DelIns Validator."""

    def validator_instance(self):
        """Return coding DNA delins instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        return CodingDNADelIns(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta
        )

    def classifier_instance(self):
        """Return the coding DNA delins classifier instance."""
        return CodingDNADelInsClassifier()

    def fixture_name(self):
        """Return the fixture name for coding DNA delins."""
        return 'coding_dna_delins'
