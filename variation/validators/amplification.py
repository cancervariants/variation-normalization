"""Module for Amplification validation"""
from typing import List, Dict, Optional

from ga4gh.vrsatile.pydantic.vrs_models import CopyChange, CopyNumberChange
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify

from variation.schemas.token_response_schema import TokenType, GeneToken
from variation.schemas.classification_response_schema import (
    Classification, ClassificationType, AmplificationClassification, Nomenclature
)
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import Validator
from variation.utils import get_priority_sequence_location


class Amplification(Validator):
    """The Insertion Validator Base class."""

    async def get_valid_invalid_results(
        self, classification: AmplificationClassification,
        transcripts: List, gene_tokens: List[GeneToken]
    ) -> List[ValidationResult]:
        # Does not require any validation
        return [ValidationResult(
            accession=None,
            classification=classification,
            is_valid=True,
            errors=[],
            gene_tokens=gene_tokens
        )]
        # validation_results = []
        # errors = []
        # cx = None


        #     gene = gene_match_token.token
        #     gene_descriptor = gene_match_token.gene_descriptor
        #     seq_loc = get_priority_sequence_location(
        #         gene_descriptor, self.seqrepo_access)
        #     if seq_loc:
        #         seq_loc_vo = models.SequenceLocation(**seq_loc)
        #         seq_loc_vo._id = ga4gh_identify((seq_loc_vo))
        #         variation = {
        #             "type": "CopyNumberChange",
        #             "subject": seq_loc_vo.as_dict(),
        #             "copy_change": CopyChange.HIGH_LEVEL_GAIN.value
        #         }
        #         variation["_id"] = ga4gh_identify(
        #             models.CopyNumberChange(**variation)
        #         )
        #         cx = CopyNumberChange(**variation).dict(by_alias=True)
        #     else:
        #         errors.append(f"No SequenceLocation found for gene: {gene}")

        # self.add_validation_result(cx, valid_variations, results, classification,
        #                             s, None, gene_tokens, errors)

    def get_gene_tokens(self, classification: Classification) -> List[GeneToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_protein_gene_symbol_tokens(classification)

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Check that classification type can be validated by validator.

        :param ClassificationType classification_type: The type of variation
        :return: `True` if classification_type matches validator's
            classification type. `False` otherwise.
        """
        return classification_type == ClassificationType.AMPLIFICATION

    def variation_name(self) -> str:
        """Return the variation name."""
        return "amplification"

    async def get_transcripts(
        self, gene_tokens: List, classification: Classification, errors: List
    ) -> List:
        """Return empty list since amplification does not require transcripts

        :param List gene_tokens: A list of gene tokens
        :param Classification classification: A classification for a list of tokens
        :param List errors: List of errors
        :return: Empty list
        """
        return []

    def get_gene_tokens(self, classification: Classification) -> List[GeneToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return [classification.gene]
