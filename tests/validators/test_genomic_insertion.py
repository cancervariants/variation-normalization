"""Module for testing Genomic Insertion Validator."""
import unittest
from variation.validators import GenomicInsertion
from variation.classifiers import GenomicInsertionClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from gene.query import QueryHandler as GeneQueryHandler


class TestGenomicInsertionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the GenomicInsertion Validator."""

    def validator_instance(self):
        """Return genomic insertion instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        tlr = Translator(data_proxy=dp)
        gene_normalizer = GeneQueryHandler()
        return GenomicInsertion(
            seqrepo_access, transcript_mappings,
            GeneSymbol(gene_normalizer),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta, dp, tlr, gene_normalizer
        )

    def classifier_instance(self):
        """Return the genomic insertion classifier instance."""
        return GenomicInsertionClassifier()

    def fixture_name(self):
        """Return the fixture name for genomic insertion."""
        return 'genomic_insertion'
