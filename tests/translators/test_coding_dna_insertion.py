"""Module for testing Coding DNA Insertion Translator."""
import unittest
from variation.classifiers import CodingDNAInsertionClassifier
from variation.translators import CodingDNAInsertion
from variation.validators import CodingDNAInsertion as CDNAD_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from gene.query import QueryHandler as GeneQueryHandler


class TestCodingDNAInsertionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Coding DNA Insertion Translator."""

    def classifier_instance(self):
        """Return coding DNA insertion instance."""
        return CodingDNAInsertionClassifier()

    def validator_instance(self):
        """Return coding DNA delins instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        tlr = Translator(data_proxy=dp)
        return CDNAD_V(
            seqrepo_access, transcript_mappings,
            GeneSymbol(GeneQueryHandler()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta, dp, tlr
        )

    def translator_instance(self):
        """Return coding DNA insertion instance."""
        return CodingDNAInsertion()

    def fixture_name(self):
        """Return the fixture name for coding DNA insertion."""
        return 'coding_dna_insertion'
