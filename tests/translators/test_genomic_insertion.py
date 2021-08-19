"""Module for testing Genomic Insertion Translator."""
import unittest
from variation.classifiers import GenomicInsertionClassifier
from variation.translators import GenomicInsertion
from variation.validators import GenomicInsertion as GD_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from gene.query import QueryHandler as GeneQueryHandler


class TestGenomicInsertionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Genomic Insertion Translator."""

    def classifier_instance(self):
        """Return genomic insertion instance."""
        return GenomicInsertionClassifier()

    def validator_instance(self):
        """Return coding DNA insertion instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        tlr = Translator(data_proxy=dp)
        return GD_V(
            seqrepo_access, transcript_mappings,
            GeneSymbol(GeneQueryHandler()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta, dp, tlr
        )

    def translator_instance(self):
        """Return genomic insertion instance."""
        return GenomicInsertion()

    def fixture_name(self):
        """Return the fixture name for genomic insertion."""
        return 'genomic_insertion'
