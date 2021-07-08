"""Module for testing Polypeptide Truncation Validator."""
import unittest
from variation.validators import PolypeptideTruncation
from variation.classifiers import PolypeptideTruncationClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript


class TestPolypeptideTruncationValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Polypeptide Truncation Validator."""

    def validator_instance(self):
        """Return Polypeptide Truncation instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        return PolypeptideTruncation(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta,
            AminoAcidCache())

    def classifier_instance(self):
        """Return the Polypeptide Truncation classifier instance."""
        return PolypeptideTruncationClassifier()

    def fixture_name(self):
        """Return the fixture name for Polypeptide Truncation."""
        return 'polypeptide_truncation'
