"""Module for providing methods for handling queries."""
from typing import Optional

from cool_seq_tool.app import CoolSeqTool
from ga4gh.vrs.extras.translator import Translator as VrsPythonTranslator
from gene.database import create_db
from gene.query import QueryHandler as GeneQueryHandler

from variation.classify import Classify
from variation.gnomad_vcf_to_protein_variation import GnomadVcfToProteinVariation
from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.normalize import Normalize
from variation.to_copy_number_variation import ToCopyNumberVariation
from variation.to_vrs import ToVRS, VRSRepresentation
from variation.tokenize import Tokenize
from variation.tokenizers import GeneSymbol
from variation.translate import Translate
from variation.validate import Validate


class QueryHandler:
    """Class for initializing handlers that make app queries."""

    def __init__(
        self,
        gene_query_handler: Optional[GeneQueryHandler] = None,
    ) -> None:
        """Initialize QueryHandler instance.
        :param gene_query_handler: Gene normalizer query handler instance. If this is
            provided, will use a current instance. If this is not provided, will create
            a new instance.
        :param uta_db_url: URL for UTA database
        """
        cool_seq_tool = CoolSeqTool()
        self.seqrepo_access = cool_seq_tool.seqrepo_access

        if not gene_query_handler:
            gene_query_handler = GeneQueryHandler(create_db())

        vrs_representation = VRSRepresentation(self.seqrepo_access)
        gene_symbol = GeneSymbol(gene_query_handler)
        tokenizer = Tokenize(gene_symbol)
        classifier = Classify()
        uta_db = cool_seq_tool.uta_db
        self.alignment_mapper = cool_seq_tool.alignment_mapper
        mane_transcript = cool_seq_tool.mane_transcript
        transcript_mappings = cool_seq_tool.transcript_mappings
        self.vrs_python_tlr = VrsPythonTranslator(data_proxy=self.seqrepo_access)
        validator = Validate(
            self.seqrepo_access, transcript_mappings, uta_db, gene_query_handler
        )
        hgvs_dup_del_mode = HGVSDupDelMode(self.seqrepo_access)
        translator = Translate(
            self.seqrepo_access,
            mane_transcript,
            uta_db,
            vrs_representation,
            hgvs_dup_del_mode,
        )
        to_vrs_params = [
            self.seqrepo_access,
            tokenizer,
            classifier,
            validator,
            translator,
        ]
        self.to_vrs_handler = ToVRS(*to_vrs_params)
        normalize_params = to_vrs_params + [
            gene_query_handler,
            transcript_mappings,
            uta_db,
        ]
        self.normalize_handler = Normalize(*normalize_params)

        mane_transcript_mappings = cool_seq_tool.mane_transcript_mappings
        to_protein_params = normalize_params + [
            mane_transcript,
            mane_transcript_mappings,
        ]
        self.gnomad_vcf_to_protein_handler = GnomadVcfToProteinVariation(
            *to_protein_params
        )
        self.to_copy_number_handler = ToCopyNumberVariation(
            *to_vrs_params + [gene_query_handler, uta_db]
        )
