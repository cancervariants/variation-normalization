"""Module for testing Silent Mutation Validator."""
import unittest
from variation.validators import SilentMutation
from variation.classifiers import SilentMutationClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import AminoAcidCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from gene.query import QueryHandler as GeneQueryHandler


class TestSilentMutationValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Silent Mutation Validator."""

    def validator_instance(self):
        """Return Silent Mutation instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        tlr = Translator(data_proxy=dp)
        return SilentMutation(
            seqrepo_access, transcript_mappings,
            GeneSymbol(GeneQueryHandler()),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta, dp, tlr, AminoAcidCache())

    def classifier_instance(self):
        """Return the Silent Mutation classifier instance."""
        return SilentMutationClassifier()

    def fixture_name(self):
        """Return the fixture name for Silent Mutation."""
        return 'silent_mutation'
