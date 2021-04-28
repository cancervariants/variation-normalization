"""Module for testing Genomic Substitution Validator."""
import unittest
from variant.validators import GenomicSubstitution
from variant.classifiers import GenomicSubstitutionClassifier
from .validator_base import ValidatorBase
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache
from variant.data_sources import TranscriptMappings, SeqRepoAccess
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestGenomicSubstitutionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Genomic Substitution Validator."""

    def validator_instance(self):
        """Return genomic substitution instance."""
        return GenomicSubstitution(SeqRepoAccess(SEQREPO_DATA_PATH),
                                   TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),  # noqa: E501
                                   GeneSymbol(GeneSymbolCache()))

    def classifier_instance(self):
        """Return the genomic substitution classifier instance."""
        return GenomicSubstitutionClassifier()

    def fixture_name(self):
        """Return the fixture name for genomic substitution."""
        return 'genomic_substitution'
