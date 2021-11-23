"""Module for testing Genomic DelIns Validator."""
import unittest
from variation.validators import GenomicDelIns
from variation.classifiers import GenomicDelInsClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from gene.query import QueryHandler as GeneQueryHandler


class TestGenomicDelInsValidator(ValidatorBase, unittest.TestCase):
    """A class to test the GenomicDelIns Validator."""

    def validator_instance(self):
        """Return genomic delins instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        tlr = Translator(data_proxy=dp)
        gene_normalizer = GeneQueryHandler()
        return GenomicDelIns(
            seqrepo_access, transcript_mappings,
            GeneSymbol(gene_normalizer),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta, dp, tlr, gene_normalizer
        )

    def classifier_instance(self):
        """Return the genomic delins classifier instance."""
        return GenomicDelInsClassifier()

    def fixture_name(self):
        """Return the fixture name for genomic delins."""
        return 'genomic_delins'
