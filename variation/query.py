"""This module provides methods for handling queries."""
from typing import Optional
from os import environ

from gene.query import QueryHandler as GeneQueryHandler
from ga4gh.vrs.extras.translator import Translator
from cool_seq_tool import TRANSCRIPT_MAPPINGS_PATH, LRG_REFSEQGENE_PATH, \
    MANE_SUMMARY_PATH, CoolSeqTool

from variation import UTA_DB_URL
from variation.tokenizers import GeneSymbol
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

    def __init__(
        self, dynamodb_url: str = "", dynamodb_region: str = "us-east-2",
        gene_query_handler: Optional[GeneQueryHandler] = None,
        transcript_file_path: str = None, refseq_file_path: str = None,
        mane_data_path: str = None, uta_db_url: str = UTA_DB_URL,
        uta_db_pwd: Optional[str] = None
    ) -> None:
        """Initialize QueryHandler instance.
        :param str dynamodb_url: URL to gene normalizer dynamodb. Only used when
            `gene_query_handler` is `None`.
        :param str dynamodb_region: AWS region for gene normalizer db. Only used when
            `gene_query_handler` is `None`.
        :param Optional[GeneQueryHandler] gene_query_handler: Gene normalizer query
            handler instance. If this is provided, will use a current instance. If this
            is not provided, will create a new instance.
        :param str transcript_file_path: Path to transcript mappings file
        :param str refseq_file_path: Path to refseq gene symbol file
        :param str mane_data_path: Path to refseq mane data file
        :param str uta_db_url: URL for UTA database
        :param Optional[str] uta_db_pwd: Password for UTA database user
        """
        if not transcript_file_path:
            transcript_file_path = environ.get("TRANSCRIPT_MAPPINGS_PATH",
                                               TRANSCRIPT_MAPPINGS_PATH)
        if not refseq_file_path:
            refseq_file_path = environ.get("LRG_REFSEQGENE_PATH", LRG_REFSEQGENE_PATH)

        if not mane_data_path:
            mane_data_path = environ.get("MANE_SUMMARY_PATH", MANE_SUMMARY_PATH)

        cool_seq_tool = CoolSeqTool(
            transcript_file_path=transcript_file_path,
            lrg_refseqgene_path=refseq_file_path,
            mane_data_path=mane_data_path,
            db_url=uta_db_url,
            db_pwd=uta_db_pwd,
            gene_query_handler=gene_query_handler,
            gene_db_url=dynamodb_url,
            gene_db_region=dynamodb_region
        )
        self._seqrepo_access = cool_seq_tool.seqrepo_access
        gene_query_handler = cool_seq_tool.gene_query_handler

        vrs_representation = VRSRepresentation(self._seqrepo_access)
        gene_symbol = GeneSymbol(gene_query_handler)
        tokenizer = Tokenize(gene_symbol)
        classifier = Classify()
        uta_db = cool_seq_tool.uta_db
        self.alignment_mapper = cool_seq_tool.alignment_mapper
        mane_transcript = cool_seq_tool.mane_transcript
        transcript_mappings = cool_seq_tool.transcript_mappings
        self._tlr = Translator(data_proxy=self._seqrepo_access)
        validator = Validate(
            self._seqrepo_access, transcript_mappings, gene_symbol, mane_transcript,
            uta_db, self._tlr, gene_query_handler, vrs_representation
        )
        translator = Translate()
        hgvs_dup_del_mode = HGVSDupDelMode(self._seqrepo_access)
        to_vrs_params = [self._seqrepo_access, tokenizer,
                         classifier, validator, translator, hgvs_dup_del_mode,
                         gene_query_handler]
        self.to_vrs_handler = ToVRS(*to_vrs_params)
        self.to_vrsatile_handler = ToVRSATILE(*to_vrs_params)
        self.normalize_handler = Normalize(*to_vrs_params + [uta_db])

        codon_table = CodonTable()
        mane_transcript_mappings = cool_seq_tool.mane_transcript_mappings
        to_protein_params = to_vrs_params + [uta_db, mane_transcript,
                                             mane_transcript_mappings, codon_table]
        self.gnomad_vcf_to_protein_handler = GnomadVcfToProteinVariation(
            *to_protein_params)
        self.to_canonical_handler = ToCanonicalVariation(
            *to_vrs_params + [self._tlr, uta_db])
        self.to_copy_number_handler = ToCopyNumberVariation(*to_vrs_params)
