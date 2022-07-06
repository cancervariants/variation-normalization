"""This module provides methods for handling queries."""
from typing import Optional
from os import environ

from gene.query import QueryHandler as GeneQueryHandler
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from uta_tools import SEQREPO_DATA_PATH, TRANSCRIPT_MAPPINGS_PATH, \
    LRG_REFSEQGENE_PATH, MANE_SUMMARY_PATH, UTATools

from variation import AMINO_ACID_PATH, UTA_DB_URL
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import AminoAcidCache
from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.to_vrs import VRSRepresentation, ToVRS
from variation.classifiers import Classify
from variation.tokenizers import Tokenize
from variation.validators import Validate
from variation.translators import Translate
from variation.data_sources import CodonTable
from variation.to_vrsatile import ToVRSATILE
from variation.gnomad_vcf_to_protein_variation import GnomadVcfToProteinVariation
from variation.normalize import Normalize
from variation.to_canonical_variation import ToCanonicalVariation
from variation.to_copy_number_variation import ToCopyNumberVariation


class QueryHandler:
    """Class for initializing handlers that make app queries."""

    def __init__(self,
                 dynamodb_url: str = "",
                 dynamodb_region: str = "us-east-2",
                 seqrepo_data_path: str = None,
                 amino_acids_file_path: str = AMINO_ACID_PATH,
                 transcript_file_path: str = None,
                 refseq_file_path: str = None,
                 mane_data_path: str = None,
                 uta_db_url: str = UTA_DB_URL,
                 uta_db_pwd: Optional[str] = None) -> None:
        """Initialize QueryHandler instance.
        :param str dynamodb_url: URL to gene-normalizer database source.
        :param str dynamodb_region: AWS default region for gene-normalizer.
        :param str seqrepo_data_path: Path to seqrepo data directory
        :param str amino_acids_file_path: Path to amino acids file
        :param str transcript_file_path: Path to transcript mappings file
        :param str refseq_file_path: Path to refseq gene symbol file
        :param str mane_data_path: Path to refseq mane data file
        :param str uta_db_url: URL for UTA database
        :param Optional[str] uta_db_pwd: Password for UTA database user
        """
        if not seqrepo_data_path:
            seqrepo_data_path = environ.get("SEQREPO_DATA_PATH", SEQREPO_DATA_PATH)
        if not transcript_file_path:
            transcript_file_path = environ.get("TRANSCRIPT_MAPPINGS_PATH",
                                               TRANSCRIPT_MAPPINGS_PATH)
        if not refseq_file_path:
            refseq_file_path = environ.get("LRG_REFSEQGENE_PATH", LRG_REFSEQGENE_PATH)
        if not mane_data_path:
            mane_data_path = environ.get("MANE_SUMMARY_PATH", MANE_SUMMARY_PATH)
        uta_tools = UTATools(seqrepo_data_path=seqrepo_data_path,
                             transcript_file_path=transcript_file_path,
                             lrg_refseqgene_path=refseq_file_path,
                             mane_data_path=mane_data_path,
                             db_url=uta_db_url, db_pwd=uta_db_pwd)
        self._seqrepo_access = uta_tools.seqrepo_access

        dp = SeqRepoDataProxy(self._seqrepo_access.seqrepo_client)

        vrs_representation = VRSRepresentation(dp, self._seqrepo_access)

        amino_acid_cache = AminoAcidCache(amino_acids_file_path=amino_acids_file_path)
        gene_normalizer = GeneQueryHandler(db_url=dynamodb_url,
                                           db_region=dynamodb_region)
        gene_symbol = GeneSymbol(gene_normalizer)
        tokenizer = Tokenize(amino_acid_cache, gene_symbol)
        classifier = Classify()
        uta_db = uta_tools.uta_db
        mane_transcript = uta_tools.mane_transcript
        transcript_mappings = uta_tools.transcript_mappings
        self._tlr = Translator(data_proxy=dp)
        validator = Validate(
            self._seqrepo_access, transcript_mappings, gene_symbol, mane_transcript,
            uta_db, self._tlr, amino_acid_cache, gene_normalizer, vrs_representation
        )
        translator = Translate()
        hgvs_dup_del_mode = HGVSDupDelMode(self._seqrepo_access)
        to_vrs_params = [self._seqrepo_access, dp, tokenizer, classifier, validator,
                         translator, hgvs_dup_del_mode]
        self.to_vrs_handler = ToVRS(*to_vrs_params)
        self.to_vrsatile_handler = ToVRSATILE(*to_vrs_params + [gene_normalizer])
        self.normalize_handler = Normalize(*to_vrs_params + [gene_normalizer, uta_db])

        codon_table = CodonTable(amino_acid_cache)
        mane_transcript_mappings = uta_tools.mane_transcript_mappings
        to_protein_params = to_vrs_params + [gene_normalizer, uta_db,
                                             mane_transcript, mane_transcript_mappings,
                                             codon_table]
        self.gnomad_vcf_to_protein_handler = GnomadVcfToProteinVariation(
            *to_protein_params)
        self.to_canonical_handler = ToCanonicalVariation(
            *to_vrs_params + [self._tlr, uta_db])
        self.to_copy_number_handler = ToCopyNumberVariation(*to_vrs_params)
