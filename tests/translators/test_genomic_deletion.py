"""Module for testing Genomic Deletion Translator."""
import unittest
from variation.classifiers import GenomicDeletionClassifier
from variation.translators import GenomicDeletion
from variation.validators import GenomicDeletion as GD_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator


class TestGenomicDeletionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Genomic Deletion Translator."""

    def classifier_instance(self):
        """Return genomic deletion instance."""
        return GenomicDeletionClassifier()

    def validator_instance(self):
        """Return coding DNA deletion instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        tlr = Translator(data_proxy=dp)
        return GD_V(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta, dp, tlr
        )

    def translator_instance(self):
        """Return genomic deletion instance."""
        return GenomicDeletion()

    def fixture_name(self):
        """Return the fixture name for genomic deletion."""
        return 'genomic_deletion'
