"""The module for Genomic Deletion Range Validation."""
import logging
from typing import List, Optional, Dict, Tuple

from ga4gh.vrsatile.pydantic.vrs_models import RelativeCopyClass
from ga4gh.vrs.extras.translator import Translator
from ga4gh.vrs import models
from gene.query import QueryHandler as GeneQueryHandler
from uta_tools.data_sources import SeqRepoAccess, TranscriptMappings, UTADatabase, \
    MANETranscript

from variation.schemas.classification_response_schema import \
    ClassificationType, Classification
from variation.schemas.app_schemas import Endpoint
from variation.schemas.token_response_schema import Token
from variation.schemas.token_response_schema import GeneMatchToken
from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.tokenizers import GeneSymbol
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum
from variation.validators.duplication_deletion_base import\
    DuplicationDeletionBase
from variation.vrs_representation import VRSRepresentation

logger = logging.getLogger("variation")
logger.setLevel(logging.DEBUG)


class GenomicDeletionRange(DuplicationDeletionBase):
    """The Genomic Deletion Range Validator class."""

    def __init__(self, seq_repo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol,
                 mane_transcript: MANETranscript,
                 uta: UTADatabase, tlr: Translator,
                 gene_normalizer: GeneQueryHandler, vrs: VRSRepresentation) -> None:
        """Initialize the Genomic Deletion Range validator.

        :param SeqRepoAccess seq_repo_access: Access to SeqRepo data
        :param TranscriptMappings transcript_mappings: Access to transcript
            mappings
        :param GeneSymbol gene_symbol: Gene symbol tokenizer
        :param MANETranscript mane_transcript: Access MANE Transcript
            information
        :param UTADatabase uta: Access to UTA queries
        :param Translator tlr: Class for translating nomenclatures to and from VRS
        :param GeneQueryHandler gene_normalizer: Access to gene-normalizer
        :param VRSRepresentation vrs: Class for creating VRS objects
        """
        super().__init__(
            seq_repo_access, transcript_mappings, gene_symbol, mane_transcript,
            uta, tlr, gene_normalizer, vrs
        )
        self.hgvs_dup_del_mode = HGVSDupDelMode(seq_repo_access)

    async def get_transcripts(self, gene_tokens: List, classification: Classification,
                              errors: List) -> Optional[List[str]]:
        """Get transcript accessions for a given classification.

        :param List gene_tokens: A list of gene tokens
        :param Classification classification: A classification for a list of
            tokens
        :param List errors: List of errors
        :return: List of transcript accessions
        """
        transcripts = await self.get_genomic_transcripts(classification, errors)
        return transcripts

    async def get_valid_invalid_results(
        self, classification_tokens: List, transcripts: List,
        classification: Classification, results: List, gene_tokens: List,
        mane_data_found: Dict, is_identifier: bool,
        hgvs_dup_del_mode: HGVSDupDelModeEnum,
        endpoint_name: Optional[Endpoint] = None,
        baseline_copies: Optional[int] = None,
        relative_copy_class: Optional[RelativeCopyClass] = None,
        do_liftover: bool = False
    ) -> None:
        """Add validation result objects to a list of results.

        :param List classification_tokens: A list of classification Tokens
        :param List transcripts: A list of transcript accessions
        :param Classification classification: A classification for a list of
            tokens
        :param List results: Stores validation result objects
        :param List gene_tokens: List of GeneMatchTokens for a classification
        :param Dict mane_data_found: MANE Transcript information found
        :param bool is_identifier: `True` if identifier is given for exact
            location. `False` otherwise.
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Must be: `default`, `absolute_cnv`,
            `relative_cnv`, `repeated_seq_expr`, `literal_seq_expr`. This parameter
            determines how to represent HGVS dup/del expressions as VRS objects.
        :param Optional[Endpoint] endpoint_name: Then name of the endpoint being used
        :param Optional[int] baseline_copies: Baseline copies number
        :param Optional[RelativeCopyClass] relative_copy_class: The relative copy class
        :param bool do_liftover: Whether or not to liftover to GRCh38 assembly
        """
        valid_alleles = list()
        for s in classification_tokens:
            for t in transcripts:
                errors = list()
                t = self.get_accession(t, classification)

                variation = await self._get_variation(
                    s, t, errors, gene_tokens, hgvs_dup_del_mode,
                    relative_copy_class=relative_copy_class,
                    baseline_copies=baseline_copies)

                if not errors and (endpoint_name == Endpoint.NORMALIZE or do_liftover):
                    await self._get_normalize_variation(
                        gene_tokens, s, t, errors, hgvs_dup_del_mode,
                        mane_data_found, relative_copy_class=relative_copy_class,
                        baseline_copies=baseline_copies)

                self.add_validation_result(
                    variation, valid_alleles, results,
                    classification, s, t, gene_tokens, errors
                )

                if is_identifier:
                    break

        if endpoint_name == Endpoint.NORMALIZE or do_liftover:
            self.add_mane_to_validation_results(
                mane_data_found, valid_alleles, results,
                classification, gene_tokens
            )

    async def _get_variation(
            self, s: Token, t: str, errors: List, gene_tokens: List,
            hgvs_dup_del_mode: HGVSDupDelModeEnum,
            relative_copy_class: Optional[RelativeCopyClass] = None,
            baseline_copies: Optional[int] = None) -> Optional[Dict]:
        """Get variation data.

        :param Token s: Classification token
        :param str t: Accession
        :param List errors: List of errors
        :param List gene_tokens: List of GeneMatchTokens for a classification
        :param HGVSDupDelMode hgvs_dup_del_mode: Mode to use for interpreting
            HGVS duplications and deletions
        :param Optional[RelativeCopyClass] relative_copy_class: The relative copy class
        :param Optional[int] baseline_copies: Baseline copies number
        :return: Dictionary containing start/end position changes and variation
        """
        variation, start, end = None, None, None
        ival, grch38 = await self._get_ival(t, s, errors, gene_tokens)

        if not errors:
            if grch38:
                t = grch38["ac"]

            allele = self.vrs.to_vrs_allele_ranges(
                t, s.coordinate_type, s.alt_type, errors, ival)

            if start is not None and end is not None:
                pos = (start, end)
            else:
                pos = None

            variation = self.hgvs_dup_del_mode.interpret_variation(
                s.alt_type, allele, errors,
                hgvs_dup_del_mode, pos=pos, relative_copy_class=relative_copy_class,
                baseline_copies=baseline_copies)
        return variation

    async def _get_normalize_variation(
            self, gene_tokens: List, s: Token, t: str, errors: List,
            hgvs_dup_del_mode: HGVSDupDelModeEnum,
            mane_data_found: Dict,
            relative_copy_class: Optional[RelativeCopyClass] = None,
            baseline_copies: Optional[int] = None) -> None:
        """Get variation that will be returned in normalize endpoint.

        :param List gene_tokens: List of gene tokens
        :param Token s: Classification token
        :param str t: Accession
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Mode to use for
            interpreting HGVS duplications and deletions
        :param Dict mane_data_found: MANE Transcript data found for given query
        :param Optional[RelativeCopyClass] relative_copy_class: The relative copy class
        :param Optional[int] baseline_copies: Baseline copies number
        """
        ival, grch38 = await self._get_ival(t, s, errors, gene_tokens, is_norm=True)
        self.add_grch38_to_mane_data(
            t, s, errors, grch38, mane_data_found, hgvs_dup_del_mode,
            ival=ival, relative_copy_class=relative_copy_class,
            baseline_copies=baseline_copies)

    async def _get_ival(
            self, t: str, s: Token, errors: List, gene_tokens: List,
            is_norm: bool = False
    ) -> Optional[Tuple[Optional[models.SequenceInterval], Optional[Dict]]]:
        """Get ival for variations with ranges.

        :param str t: Accession
        :param Token t: Classification token
        :param List errors: List of errors
        :param List gene_tokens: LIst of gene tokens
        :param bool is_norm: `True` if normalize endpoint is being used.
            `False` otherwise.
        :return: Sequence Interval and GRCh38 data if normalize endpoint
            is being used
        """
        ival, grch38, start1, start2, end1, end2 = None, None, None, None, None, None  # noqa: E501
        # (#_#)_(#_#)
        if is_norm:
            t, start1, start2, end1, end2, grch38 = await self.get_grch38_pos_ac(
                t, s.start_pos1_del, s.start_pos2_del, pos3=s.end_pos1_del,
                pos4=s.end_pos2_del
            )
        else:
            start1 = s.start_pos1_del
            start2 = s.start_pos2_del
            end1 = s.end_pos1_del
            end2 = s.end_pos2_del

        gene = gene_tokens[0].token if gene_tokens else None
        await self.validate_gene_or_accession_pos(
            t, [start1, start2, end1, end2], errors, gene=gene)

        if not errors and start1 and start2 and end1 and end2:
            ival = self.vrs.get_ival_certain_range(
                start1, start2, end1, end2
            )

        return ival, grch38

    def get_gene_tokens(
            self, classification: Classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_gene_symbol_tokens(classification)

    def variation_name(self) -> str:
        """Return the variation name."""
        return "genomic deletion range"

    def is_token_instance(self, t: Token) -> bool:
        """Check that token is Genomic Deletion Range.

        :param Token t: Classification token
        """
        return t.token_type == "GenomicDeletionRange"

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is
        Genomic Deletion Range.

        :param ClassificationType classification_type: Classification type
        :return: `True` if classification type matches, `False` otherwise
        """
        return classification_type == \
            ClassificationType.GENOMIC_DELETION_RANGE
