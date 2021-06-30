"""Module for testing Genomic Substitution Validator."""
import unittest
from variation.validators import GenomicSubstitution
from variation.classifiers import GenomicSubstitutionClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript


class TestGenomicSubstitutionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Genomic Substitution Validator."""

    def validator_instance(self):
        """Return genomic substitution instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        return GenomicSubstitution(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), UTA())
        )

    def classifier_instance(self):
        """Return the genomic substitution classifier instance."""
        return GenomicSubstitutionClassifier()

    def fixture_name(self):
        """Return the fixture name for genomic substitution."""
        return 'genomic_substitution'
