"""Module for testing amino acid deletion Translator."""
import unittest
from variation.classifiers import AminoAcidDeletionClassifier
from variation.translators import AminoAcidDeletion
from variation.validators import AminoAcidDeletion as AAD_V
from .translator_base import TranslatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator


class TestAminoAcidDeletionTranslator(TranslatorBase, unittest.TestCase):
    """A class to test the amino acid deletion Translator."""

    def classifier_instance(self):
        """Return amino acid deletion instance."""
        return AminoAcidDeletionClassifier()

    def validator_instance(self):
        """Return amino acid deletion instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        tlr = Translator(data_proxy=dp)
        return AAD_V(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta, dp, tlr, AminoAcidCache())

    def translator_instance(self):
        """Return amino acid deletion instance."""
        return AminoAcidDeletion()

    def fixture_name(self):
        """Return the fixture name for amino acid deletion."""
        return 'amino_acid_deletion'
