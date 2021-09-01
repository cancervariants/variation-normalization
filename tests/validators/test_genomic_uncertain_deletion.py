"""Module for testing Genomic Uncertain Deletion Validator."""
import unittest
from variation.validators import GenomicUncertainDeletion
from variation.classifiers import GenomicUncertainDeletionClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from gene.query import QueryHandler as GeneQueryHandler


class TestGenomicUncertainDeletionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Genomic Uncertain Deletion Validator."""

    def validator_instance(self):
        """Return genomic uncertain deletion instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        tlr = Translator(data_proxy=dp)
        return GenomicUncertainDeletion(
            seqrepo_access, transcript_mappings,
            GeneSymbol(GeneQueryHandler()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta, dp, tlr
        )

    def classifier_instance(self):
        """Return the genomic uncertain deletion classifier instance."""
        return GenomicUncertainDeletionClassifier()

    def fixture_name(self):
        """Return the fixture name for genomic uncertain deletion."""
        return 'genomic_uncertain_deletion'
