"""Module for testing Coding DNA Silent Mutation Translator."""
import unittest
from variation.classifiers import CodingDNASilentMutationClassifier
from variation.translators import CodingDNASilentMutation
from variation.validators import CodingDNASilentMutation as CDNASM_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from gene.query import QueryHandler as GeneQueryHandler


class TestCodingDNASilentMutationTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Coding DNA Silent Mutation Translator."""

    def classifier_instance(self):
        """Return coding DNA silent mutation instance."""
        return CodingDNASilentMutationClassifier()

    def validator_instance(self):
        """Return coding DNA silent mutation instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        tlr = Translator(data_proxy=dp)
        return CDNASM_V(
            seqrepo_access, transcript_mappings,
            GeneSymbol(GeneQueryHandler()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta, dp, tlr
        )

    def translator_instance(self):
        """Return coding DNA silent mutation instance."""
        return CodingDNASilentMutation()

    def fixture_name(self):
        """Return the fixture name for coding DNA silent mutation."""
        return 'coding_dna_silent_mutation'
