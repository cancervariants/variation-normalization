"""Module for testing Genomic Silent Mutation Validator."""
import unittest
from variation.validators import GenomicSilentMutation
from variation.classifiers import GenomicSilentMutationClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import GeneSymbolCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator


class TestGenomicSilentMutationValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Genomic Silent Mutation Validator."""

    def validator_instance(self):
        """Return genomic silent mutation instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        tlr = Translator(data_proxy=dp)
        return GenomicSilentMutation(
            seqrepo_access, transcript_mappings, GeneSymbol(GeneSymbolCache()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta, dp, tlr
        )

    def classifier_instance(self):
        """Return the genomic silent mutation classifier instance."""
        return GenomicSilentMutationClassifier()

    def fixture_name(self):
        """Return the fixture name for genomic silent mutation."""
        return 'genomic_silent_mutation'
