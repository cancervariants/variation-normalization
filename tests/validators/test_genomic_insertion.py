"""Module for testing Genomic Insertion Validator."""
import unittest
from variation.validators import GenomicInsertion
from variation.classifiers import GenomicInsertionClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript


class TestGenomicInsertionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the GenomicInsertion Validator."""

    def validator_instance(self):
        """Return genomic insertion instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        return GenomicInsertion(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta
        )

    def classifier_instance(self):
        """Return the genomic insertion classifier instance."""
        return GenomicInsertionClassifier()

    def fixture_name(self):
        """Return the fixture name for genomic insertion."""
        return 'genomic_insertion'
