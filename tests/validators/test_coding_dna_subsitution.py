"""Module for testing Coding DNA Substitution Validator."""
import unittest
from variation.validators import CodingDNASubstitution
from variation.classifiers import CodingDNASubstitutionClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from gene.query import QueryHandler as GeneQueryHandler


class TestCodingDNASubstitutionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Coding DNA Substitution Validator."""

    def validator_instance(self):
        """Return coding DNA substitution instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        tlr = Translator(data_proxy=dp)
        return CodingDNASubstitution(
            seqrepo_access, transcript_mappings,
            GeneSymbol(GeneQueryHandler()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta, dp, tlr
        )

    def classifier_instance(self):
        """Return the coding DNA substitution classifier instance."""
        return CodingDNASubstitutionClassifier()

    def fixture_name(self):
        """Return the fixture name for coding DNA substitution."""
        return 'coding_dna_substitution'
