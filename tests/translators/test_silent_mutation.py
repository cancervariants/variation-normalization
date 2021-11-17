"""Module for testing Silent Mutation Translator."""
import unittest
from variation.classifiers import SilentMutationClassifier
from variation.translators import SilentMutation
from variation.validators import SilentMutation as SM_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import AminoAcidCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from gene.query import QueryHandler as GeneQueryHandler


class TestSilentMutationTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the silent mutation Translator."""

    def classifier_instance(self):
        """Return silent mutation instance."""
        return SilentMutationClassifier()

    def validator_instance(self):
        """Return silent mutation instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        tlr = Translator(data_proxy=dp)
        gene_normalizer = GeneQueryHandler()
        return SM_V(
            seqrepo_access, transcript_mappings,
            GeneSymbol(gene_normalizer),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta, dp, tlr, gene_normalizer, AminoAcidCache())

    def translator_instance(self):
        """Return silent mutation instance."""
        return SilentMutation()

    def fixture_name(self):
        """Return the fixture name for silent mutation."""
        return 'silent_mutation'
