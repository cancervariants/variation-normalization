"""Module for testing Amino Acid DelIns Validator."""
import unittest
from variation.validators import AminoAcidDelIns
from variation.classifiers import AminoAcidDelInsClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript


class TestAminoAcidDelInsValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Protein DelIns Validator."""

    def validator_instance(self):
        """Return amino acid delins instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        return AminoAcidDelIns(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), UTA()),
            AminoAcidCache())

    def classifier_instance(self):
        """Return the amino acid delins classifier instance."""
        return AminoAcidDelInsClassifier()

    def fixture_name(self):
        """Return the fixture name for amino acid delins."""
        return 'amino_acid_delins'
