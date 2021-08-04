"""Module for testing Genomic Silent Mutation Translator."""
import unittest
from variation.classifiers import GenomicSilentMutationClassifier
from variation.translators import GenomicSilentMutation
from variation.validators import GenomicSilentMutation as GENOMICSM_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator


class TestGenomicSilentMutationTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Genomic Silent Mutation Translator."""

    def classifier_instance(self):
        """Return genomic silent mutation instance."""
        return GenomicSilentMutationClassifier()

    def validator_instance(self):
        """Return genomic silent mutation instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        tlr = Translator(data_proxy=dp)
        return GENOMICSM_V(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta, dp, tlr
        )

    def translator_instance(self):
        """Return genomic silent mutation instance."""
        return GenomicSilentMutation()

    def fixture_name(self):
        """Return the fixture name for genomic silent mutation."""
        return 'genomic_silent_mutation'
