"""Module for testing Polypeptide Truncation Validator."""
import unittest
from variation.validators import PolypeptideTruncation
from variation.classifiers import PolypeptideTruncationClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import AminoAcidCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from gene.query import QueryHandler as GeneQueryHandler


class TestPolypeptideTruncationValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Polypeptide Truncation Validator."""

    def validator_instance(self):
        """Return Polypeptide Truncation instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        tlr = Translator(data_proxy=dp)
        return PolypeptideTruncation(
            seqrepo_access, transcript_mappings,
            GeneSymbol(GeneQueryHandler()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta, dp, tlr, AminoAcidCache())

    def classifier_instance(self):
        """Return the Polypeptide Truncation classifier instance."""
        return PolypeptideTruncationClassifier()

    def fixture_name(self):
        """Return the fixture name for Polypeptide Truncation."""
        return 'polypeptide_truncation'
