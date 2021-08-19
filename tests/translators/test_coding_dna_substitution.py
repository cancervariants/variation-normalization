"""Module for testing Coding DNA Substitution Translator."""
import unittest
from variation.classifiers import CodingDNASubstitutionClassifier
from variation.translators import CodingDNASubstitution
from variation.validators import CodingDNASubstitution as CDNASUB_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from gene.query import QueryHandler as GeneQueryHandler


class TestCodingDNASubstitutionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the Coding DNA Substitution Translator."""

    def classifier_instance(self):
        """Return coding DNA substitution instance."""
        return CodingDNASubstitutionClassifier()

    def validator_instance(self):
        """Return coding DNA substitution instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        tlr = Translator(data_proxy=dp)
        return CDNASUB_V(
            seqrepo_access, transcript_mappings,
            GeneSymbol(GeneQueryHandler()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta, dp, tlr
        )

    def translator_instance(self):
        """Return coding DNA substitution instance."""
        return CodingDNASubstitution()

    def fixture_name(self):
        """Return the fixture name for coding DNA substitution."""
        return 'coding_dna_substitution'
