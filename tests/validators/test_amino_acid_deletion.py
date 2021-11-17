"""Module for testing Amino Acid Deletion Validator."""
import unittest
from variation.validators import AminoAcidDeletion
from variation.classifiers import AminoAcidDeletionClassifier
from .validator_base import ValidatorBase
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import AminoAcidCache
from variation.data_sources import TranscriptMappings, SeqRepoAccess, \
    MANETranscriptMappings, UTA
from variation.mane_transcript import MANETranscript
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from gene.query import QueryHandler as GeneQueryHandler


class TestAminoAcidDeletionValidator(ValidatorBase, unittest.TestCase):
    """A class to test the Amino Acid Deletion Validator."""

    def validator_instance(self):
        """Return amino acid deletion instance."""
        seqrepo_access = SeqRepoAccess()
        transcript_mappings = TranscriptMappings()
        uta = UTA()
        dp = SeqRepoDataProxy(seqrepo_access.seq_repo_client)
        tlr = Translator(data_proxy=dp)
        gene_normalizer = GeneQueryHandler()
        return AminoAcidDeletion(
            seqrepo_access, transcript_mappings,
            GeneSymbol(gene_normalizer),
            MANETranscript(seqrepo_access, transcript_mappings,
                           MANETranscriptMappings(), uta),
            uta, dp, tlr, gene_normalizer, AminoAcidCache())

    def classifier_instance(self):
        """Return the amino acid deletion classifier instance."""
        return AminoAcidDeletionClassifier()

    def fixture_name(self):
        """Return the fixture name for amino acid deletion."""
        return 'amino_acid_deletion'
