"""Module for testing Genomic Insertion Validator."""
import unittest
from variant.validators import GenomicInsertion
from variant.classifiers import GenomicInsertionClassifier
from .validator_base import ValidatorBase
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache
from variant.data_sources import TranscriptMappings, SeqRepoAccess
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestGenomicInsertionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the GenomicInsertion Validator."""

    def validator_instance(self):
        """Return genomic insertion instance."""
        return GenomicInsertion(SeqRepoAccess(SEQREPO_DATA_PATH),
                                TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),
                                GeneSymbol(GeneSymbolCache()))

    def classifier_instance(self):
        """Return the genomic insertion classifier instance."""
        return GenomicInsertionClassifier()

    def fixture_name(self):
        """Return the fixture name for genomic insertion."""
        return 'genomic_insertion'
