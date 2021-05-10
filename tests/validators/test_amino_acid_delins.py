"""Module for testing Amino Acid DelIns Validator."""
import unittest
from variant.validators import AminoAcidDelIns
from variant.classifiers import AminoAcidDelInsClassifier
from .validator_base import ValidatorBase
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from variant.data_sources import TranscriptMappings, SeqRepoAccess
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestAminoAcidDelInsValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Protein DelIns Validator."""

    def validator_instance(self):
        """Return amino acid delins instance."""
        return AminoAcidDelIns(SeqRepoAccess(SEQREPO_DATA_PATH),
                               TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),
                               GeneSymbol(GeneSymbolCache()),
                               AminoAcidCache())

    def classifier_instance(self):
        """Return the amino acid delins classifier instance."""
        return AminoAcidDelInsClassifier()

    def fixture_name(self):
        """Return the fixture name for amino acid delins."""
        return 'amino_acid_delins'
