"""Module for testing Genomic DelIns Validator."""
import unittest
from variant.validators import GenomicDelIns
from variant.classifiers import GenomicDelInsClassifier
from .validator_base import ValidatorBase
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache
from variant.data_sources import TranscriptMappings, SeqRepoAccess
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestGenomicDelInsValidator(ValidatorBase, unittest.TestCase):
    """A class to test the GenomicDelIns Validator."""

    def validator_instance(self):
        """Return genomic delins instance."""
        return GenomicDelIns(SeqRepoAccess(SEQREPO_DATA_PATH),
                             TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),  # noqa: E501
                             GeneSymbol(GeneSymbolCache()))

    def classifier_instance(self):
        """Return the genomic delins classifier instance."""
        return GenomicDelInsClassifier()

    def fixture_name(self):
        """Return the fixture name for genomic delins."""
        return 'genomic_delins'
