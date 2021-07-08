"""Module for testing Amino Acid Insertion Validator."""
import unittest
from variation.validators import AminoAcidInsertion
from variation.classifiers import AminoAcidInsertionClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript


class TestAminoAcidInsertionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Amino Acid Insertion Validator."""

    def validator_instance(self):
        """Return amino acid insertion instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        return AminoAcidInsertion(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta,
            AminoAcidCache())

    def classifier_instance(self):
        """Return the amino acid insertion classifier instance."""
        return AminoAcidInsertionClassifier()

    def fixture_name(self):
        """Return the fixture name for amino acid insertion."""
        return 'amino_acid_insertion'
