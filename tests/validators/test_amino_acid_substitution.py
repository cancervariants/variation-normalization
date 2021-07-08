"""Module for testing Amino Acid Substitution Validator."""
import unittest
from variation.validators import AminoAcidSubstitution
from variation.classifiers import AminoAcidSubstitutionClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript


class TestAminoAcidSubstitutionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Protein Substitution Validator."""

    def validator_instance(self):
        """Return amino acid substitution instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        return AminoAcidSubstitution(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta,
            AminoAcidCache())

    def classifier_instance(self):
        """Return the protein substitution classifier instance."""
        return AminoAcidSubstitutionClassifier()

    def fixture_name(self):
        """Return the fixture name for protein substitution."""
        return 'amino_acid_substitution'
