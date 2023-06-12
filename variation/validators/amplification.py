"""Module for Amplification validation"""
from typing import List

from variation.schemas.token_response_schema import GeneToken
from variation.schemas.classification_response_schema import (
    Classification, ClassificationType, AmplificationClassification
)
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import Validator


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
