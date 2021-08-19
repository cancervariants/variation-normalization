"""Module for testing Coding DNA Deletion Validator."""
import unittest
from variation.validators import CodingDNADeletion
from variation.classifiers import CodingDNADeletionClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from gene.query import QueryHandler as GeneQueryHandler


class TestCodingDNADeletionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the CodingDNADeletion Validator."""

    def validator_instance(self):
        """Return coding dna deletion instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        tlr = Translator(data_proxy=dp)
        return CodingDNADeletion(
            seqrepo_access, transcript_mappings,
            GeneSymbol(GeneQueryHandler()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta, dp, tlr
        )

    def classifier_instance(self):
        """Return the coding dna deletion classifier instance."""
        return CodingDNADeletionClassifier()

    def fixture_name(self):
        """Return the fixture name for coding dna deletion."""
        return 'coding_dna_deletion'
