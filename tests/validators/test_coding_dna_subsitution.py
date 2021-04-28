"""Module for testing Coding DNA Substitution Validator."""
import unittest
from variant.validators import CodingDNASubstitution
from variant.classifiers import CodingDNASubstitutionClassifier
from .validator_base import ValidatorBase
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache
from variant.data_sources import TranscriptMappings, SeqRepoAccess
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestCodingDNASubstitutionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Coding DNA Substitution Validator."""

    def validator_instance(self):
        """Return coding DNA substitution instance."""
        return CodingDNASubstitution(SeqRepoAccess(SEQREPO_DATA_PATH),
                                     TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),  # noqa: E501
                                     GeneSymbol(GeneSymbolCache()))

    def classifier_instance(self):
        """Return the coding DNA substitution classifier instance."""
        return CodingDNASubstitutionClassifier()

    def fixture_name(self):
        """Return the fixture name for coding DNA substitution."""
        return 'coding_dna_substitution'
