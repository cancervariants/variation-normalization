"""Module for testing polypeptide truncation Translator."""
import unittest
from variation.classifiers import PolypeptideTruncationClassifier
from variation.translators import PolypeptideTruncation
from variation.validators import PolypeptideTruncation as PT_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator


class TestPolypeptideTruncationTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the polypeptide truncation Translator."""

    def classifier_instance(self):
        """Return polypeptide truncation instance."""
        return PolypeptideTruncationClassifier()

    def validator_instance(self):
        """Return polypeptide truncation instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        tlr = Translator(data_proxy=dp)
        return PT_V(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta, dp, tlr, AminoAcidCache())

    def translator_instance(self):
        """Return polypeptide truncation instance."""
        return PolypeptideTruncation()

    def fixture_name(self):
        """Return the fixture name for polypeptide truncation."""
        return 'polypeptide_truncation'
