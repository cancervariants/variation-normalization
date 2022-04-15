"""The module for Genomic Silent Mutation Validation."""
import logging
from typing import Optional, List, Dict

from ga4gh.vrsatile.pydantic.vrs_models import RelativeCopyClass

from variation.schemas.app_schemas import Endpoint
from variation.schemas.classification_response_schema import ClassificationType, \
    Classification
from variation.schemas.token_response_schema import Token
from variation.schemas.normalize_response_schema import HGVSDupDelMode as HGVSDupDelModeEnum   # noqa: E501
from .single_nucleotide_variation_base import SingleNucleotideVariationBase


logger = logging.getLogger("variation")
logger.setLevel(logging.DEBUG)


class GenomicSilentMutation(SingleNucleotideVariationBase):
    """The Genomic Silent Mutation Validator class."""

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
        await self.silent_mutation_valid_invalid_results(
            classification_tokens, transcripts, classification, results,
            gene_tokens, endpoint_name, mane_data_found, is_identifier,
            do_liftover=do_liftover
        )

    def get_gene_tokens(self, classification: Classification) -> List:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_gene_symbol_tokens(classification)

    def variation_name(self) -> str:
        """Return the variation name."""
        return "genomic silent mutation"

    def is_token_instance(self, t: Token) -> bool:
        """Check that token is Genomic Silent Mutation."""
        return t.token_type == "GenomicSilentMutation"

    def validates_classification_type(
            self, classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is genomic substitution."""
        return classification_type == ClassificationType.GENOMIC_SILENT_MUTATION
