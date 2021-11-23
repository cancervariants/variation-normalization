"""Module for testing Genomic Duplication Validator."""
import unittest
from variation.validators import GenomicDuplication
from variation.classifiers import GenomicDuplicationClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from gene.query import QueryHandler as GeneQueryHandler


class TestGenomicDuplicationValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Genomic Duplication Validator."""

    def validator_instance(self):
        """Return genomic duplication instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        tlr = Translator(data_proxy=dp)
        gene_normalizer = GeneQueryHandler()
        return GenomicDuplication(
            seqrepo_access, transcript_mappings,
            GeneSymbol(gene_normalizer),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta, dp, tlr, gene_normalizer
        )

    def classifier_instance(self):
        """Return the genomic duplication classifier instance."""
        return GenomicDuplicationClassifier()

    def fixture_name(self):
        """Return the fixture name for genomic duplication."""
        return 'genomic_duplication'
