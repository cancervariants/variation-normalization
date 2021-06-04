"""Module for testing Amino Acid Insertion Validator."""
import unittest
from variation.validators import AminoAcidInsertion
from variation.classifiers import AminoAcidInsertionClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess
from variation import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestAminoAcidInsertionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Amino Acid Insertion Validator."""

    def validator_instance(self):
        """Return amino acid insertion instance."""
        return AminoAcidInsertion(SeqRepoAccess(SEQREPO_DATA_PATH),
                                  TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),
                                  GeneSymbol(GeneSymbolCache()),
                                  AminoAcidCache())

    def classifier_instance(self):
        """Return the amino acid insertion classifier instance."""
        return AminoAcidInsertionClassifier()

    def fixture_name(self):
        """Return the fixture name for amino acid insertion."""
        return 'amino_acid_insertion'
