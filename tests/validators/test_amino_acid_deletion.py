"""Module for testing Amino Acid Deletion Validator."""
import unittest
from variant.validators import AminoAcidDeletion
from variant.classifiers import AminoAcidDeletionClassifier
from .validator_base import ValidatorBase
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from variant.data_sources import TranscriptMappings, SeqRepoAccess
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestAminoAcidDeletionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Amino Acid Deletion Validator."""

    def validator_instance(self):
        """Return amino acid deletion instance."""
        return AminoAcidDeletion(SeqRepoAccess(SEQREPO_DATA_PATH),
                                 TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),
                                 GeneSymbol(GeneSymbolCache()),
                                 AminoAcidCache())

    def classifier_instance(self):
        """Return the amino acid deletion classifier instance."""
        return AminoAcidDeletionClassifier()

    def fixture_name(self):
        """Return the fixture name for amino acid deletion."""
        return 'amino_acid_deletion'
