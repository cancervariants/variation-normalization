"""The module for Genomic Deletion Validation."""
import logging
from typing import List, Optional, Dict

from ga4gh.vrsatile.pydantic.vrs_models import RelativeCopyClass

from variation.schemas.app_schemas import Endpoint
from variation.validators.duplication_deletion_base import\
    DuplicationDeletionBase
from variation.schemas.classification_response_schema import \
    ClassificationType, Classification
from variation.schemas.token_response_schema import Token, SequenceOntology
from variation.schemas.token_response_schema import GeneMatchToken
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum


logger = logging.getLogger("variation")
logger.setLevel(logging.DEBUG)


class GenomicDeletion(DuplicationDeletionBase):
    """The Genomic Deletion Validator class."""

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
                start, end = s.start_pos_del, None
                if s.end_pos_del is None:
                    end = start
                else:
                    end = s.end_pos_del

                # Validate pos
                if gene_tokens:
                    await self.validate_gene_or_accession_pos(
                        t, [start, end], errors, gene=gene_tokens[0].token)
                else:
                    self.check_reference_sequence(t, s, errors)

                if not errors:
                    allele = self.vrs.to_vrs_allele(
                        t, s.start_pos_del, s.end_pos_del,
                        s.coordinate_type, s.alt_type, errors)

                    variation = self.hgvs_dup_del_mode.interpret_variation(
                        s.alt_type, allele, errors, hgvs_dup_del_mode,
                        pos=(start, end), relative_copy_class=relative_copy_class,
                        baseline_copies=baseline_copies
                    )
                else:
                    variation = None

                if not errors and (endpoint_name == Endpoint.NORMALIZE or do_liftover):
                    await self._get_normalize_variation(
                        gene_tokens, s, t, errors, hgvs_dup_del_mode,
                        mane_data_found, start, end,
                        relative_copy_class=relative_copy_class,
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

    async def _get_normalize_variation(
            self, gene_tokens: List, s: Token, t: str, errors: List,
            hgvs_dup_del_mode: HGVSDupDelModeEnum, mane_data_found: Dict,
            start: int, end: int,
            relative_copy_class: Optional[RelativeCopyClass] = None,
            baseline_copies: Optional[int] = None) -> None:
        """Get variation that will be returned in normalize endpoint.

        :param List gene_tokens: List of gene tokens
        :param Token s: Classification token
        :param str t: Accession
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Mode to use for
            interpreting HGVS duplications and deletions
        :param Dict mane_data_found: MANE Transcript data found for given query
        :param int start: Start pos change
        :param int end: End pos change
        :param Optional[RelativeCopyClass] relative_copy_class: The relative copy class
        :param Optional[int] baseline_copies: Baseline copies number
        """
        if gene_tokens:
            await self.add_normalized_genomic_dup_del(
                s, t, s.start_pos_del, s.end_pos_del, gene_tokens[0].token,
                SequenceOntology.DELETION, errors, hgvs_dup_del_mode,
                mane_data_found, relative_copy_class=relative_copy_class,
                baseline_copies=baseline_copies)
        else:
            # No gene provided, then use GRCh38 assembly
            if not self._is_grch38_assembly(t):
                grch38 = await self.mane_transcript.g_to_grch38(t, start, end)
            else:
                grch38 = dict(ac=t, pos=(start, end))

            if grch38:
                await self.validate_gene_or_accession_pos(
                    grch38["ac"], [grch38["pos"][0], grch38["pos"][1]], errors)
                self.add_grch38_to_mane_data(
                    t, s, errors, grch38, mane_data_found, hgvs_dup_del_mode,
                    use_vrs_allele_range=False, relative_copy_class=relative_copy_class,
                    baseline_copies=baseline_copies)

    def get_gene_tokens(
            self, classification: Classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_gene_symbol_tokens(classification)

    def variation_name(self) -> str:
        """Return the variation name."""
        return "genomic deletion"

    def is_token_instance(self, t: Token) -> bool:
        """Check that token is Genomic Deletion.

        :param Token t: Classification token
        """
        return t.token_type == "GenomicDeletion"

    def validates_classification_type(
            self, classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is
        Genomic DelIns.

        :param ClassificationType classification_type: Classification type
        :return: `True` if classification type matches, `False` otherwise
        """
        return classification_type == ClassificationType.GENOMIC_DELETION
