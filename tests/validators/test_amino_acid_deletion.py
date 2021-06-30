"""Module for testing Amino Acid Deletion Validator."""
import unittest
from variation.validators import AminoAcidDeletion
from variation.classifiers import AminoAcidDeletionClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript


class TestAminoAcidDeletionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Amino Acid Deletion Validator."""

    def validator_instance(self):
        """Return amino acid deletion instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        return AminoAcidDeletion(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), UTA()),
            AminoAcidCache())

    def classifier_instance(self):
        """Return the amino acid deletion classifier instance."""
        return AminoAcidDeletionClassifier()

    def fixture_name(self):
        """Return the fixture name for amino acid deletion."""
        return 'amino_acid_deletion'
