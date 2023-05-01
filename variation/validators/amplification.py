"""Module for Amplification validation"""
from typing import List, Dict, Optional

from ga4gh.vrsatile.pydantic.vrs_models import CopyChange, CopyNumberChange
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify

from variation.schemas.token_response_schema import GeneMatchToken, TokenType
from variation.schemas.classification_response_schema import Classification,\
    ClassificationType
from variation.schemas.app_schemas import Endpoint
from variation.schemas.token_response_schema import Token
from variation.validators.validator import Validator
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum
from variation.utils import get_priority_sequence_location


class Amplification(Validator):
    """The Insertion Validator Base class."""

    async def get_valid_invalid_results(
        self, classification_tokens: List, transcripts: List,
        classification: Classification, results: List, gene_tokens: List,
        mane_data_found: Dict, is_identifier: bool,
        hgvs_dup_del_mode: HGVSDupDelModeEnum,
        endpoint_name: Optional[Endpoint] = None,
        baseline_copies: Optional[int] = None,
        copy_change: Optional[CopyChange] = None,
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
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Must be: `default`,
            `copy_number_count`, `copy_number_change`, `repeated_seq_expr`,
            `literal_seq_expr`. This parameter determines how to represent HGVS dup/del
            expressions as VRS objects.
        :param Optional[Endpoint] endpoint_name: Then name of the endpoint being used
        :param Optional[int] baseline_copies: Baseline copies number
        :param Optional[CopyChange] copy_change: The copy change
        :param bool do_liftover: Whether or not to liftover to GRCh38 assembly
        """
        valid_variations = list()
        gene_match_tokens = [token for token in classification.all_tokens
                             if token.token_type == "GeneSymbol"]
        for s in classification_tokens:
            errors = list()
            cx = None

            if gene_match_tokens:
                gene_match_token = gene_match_tokens[0]
                gene = gene_match_token.token
                gene_descriptor = gene_match_token.gene_descriptor
                seq_loc = get_priority_sequence_location(
                    gene_descriptor, self.seqrepo_access)
                if seq_loc:
                    seq_loc_vo = models.SequenceLocation(**seq_loc)
                    seq_loc_vo._id = ga4gh_identify((seq_loc_vo))
                    variation = {
                        "type": "CopyNumberChange",
                        "subject": seq_loc_vo.as_dict(),
                        "copy_change": CopyChange.HIGH_LEVEL_GAIN.value
                    }
                    variation["_id"] = ga4gh_identify(
                        models.CopyNumberChange(**variation)
                    )
                    cx = CopyNumberChange(**variation).dict(by_alias=True)
                else:
                    errors.append(f"No SequenceLocation found for gene: {gene}")
            else:
                errors.append("No gene_tokens found")

            self.add_validation_result(cx, valid_variations, results, classification,
                                       s, None, gene_tokens, errors)

    def get_gene_tokens(self, classification: Classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_protein_gene_symbol_tokens(classification)

    async def get_transcripts(self, gene_tokens: List, classification: Classification,
                              errors: List) -> Optional[List[str]]:
        """Get transcript accessions for a given classification.

        :param List gene_tokens: A list of gene tokens
        :param Classification classification: A classification for a list of tokens
        :param List errors: List of errors
        :return: List of transcript accessions
        """
        return []

    def variation_name(self) -> str:
        """Return the variation name."""
        return "amplification"

    def is_token_instance(self, t: Token) -> bool:
        """Check that token is Amplification.

        :param Token t: Classification token
        :return: `True` if amplification. `False` otherwise
        """
        return t.token_type == TokenType.AMPLIFICATION

    def validates_classification_type(self,
                                      classification_type: ClassificationType) -> bool:
        """Check that classification type can be validated by validator.

        :param ClassificationType classification_type: The type of variation
        :return: `True` if classification_type matches validator's
            classification type. `False` otherwise.
        """
        return classification_type == ClassificationType.AMPLIFICATION
