"""Module for testing Genomic Substitution Translator."""
import unittest
from variation.classifiers import GenomicSubstitutionClassifier
from variation.translators import GenomicSubstitution
from variation.validators import GenomicSubstitution as GSUB_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from gene.query import QueryHandler as GeneQueryHandler


class TestGenomicSubstitutionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Genomic Substitution Translator."""

    def classifier_instance(self):
        """Return genomic substitution instance."""
        return GenomicSubstitutionClassifier()

    def validator_instance(self):
        """Return genomic substitution instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        tlr = Translator(data_proxy=dp)
        return GSUB_V(
            seqrepo_access, transcript_mappings,
            GeneSymbol(GeneQueryHandler()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta, dp, tlr
        )

    def translator_instance(self):
        """Return genomic substitution instance."""
        return GenomicSubstitution()

    def fixture_name(self):
        """Return the fixture name for genomic substitution."""
        return 'genomic_substitution'
