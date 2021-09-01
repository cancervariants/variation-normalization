"""Module for testing Genomic Uncertain Deletion Translator."""
import unittest
from variation.classifiers import GenomicUncertainDeletionClassifier
from variation.translators import GenomicUncertainDeletion
from variation.validators import GenomicUncertainDeletion as GUD_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from gene.query import QueryHandler as GeneQueryHandler


class TestGenomicUncertainDeletionTranslator(TranslatorBase,
                                             unittest.TestCase):
    """A class to test the Genomic Uncertain Deletion Translator."""

    def classifier_instance(self):
        """Return genomic uncertain deletion instance."""
        return GenomicUncertainDeletionClassifier()

    def validator_instance(self):
        """Return genomic uncertain deletion instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        tlr = Translator(data_proxy=dp)
        return GUD_V(
            seqrepo_access, transcript_mappings,
            GeneSymbol(GeneQueryHandler()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta, dp, tlr
        )

    def translator_instance(self):
        """Return genomic uncertain deletion instance."""
        return GenomicUncertainDeletion()

    def fixture_name(self):
        """Return the fixture name for genomic uncertain deletion."""
        return 'genomic_uncertain_deletion'
