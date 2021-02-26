"""Module for testing Amino Acid Substitution Validator."""
import unittest
from variant.validators import AminoAcidSubstitution
from variant.classifiers import AminoAcidSubstitutionClassifier
from .validator_base import ValidatorBase
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from variant.data_sources import TranscriptMappings, SeqRepoAccess
from variant import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH


class TestAminoAcidSubstitutionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Protein Substitution Validator."""

    def validator_instance(self):
        """Return amino acid substitution instance."""
        return AminoAcidSubstitution(SeqRepoAccess(SEQREPO_DATA_PATH),
                                     TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH),  # noqa: E501
                                     GeneSymbol(GeneSymbolCache()),
                                     AminoAcidCache())

    def classifier_instance(self):
        """Return the protein substitution classifier instance."""
        return AminoAcidSubstitutionClassifier()

    def fixture_name(self):
        """Return the fixture name for protein substitution."""
        return 'amino_acid_substitution'
