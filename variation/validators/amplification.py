"""Module for Amplification validation"""
from typing import List

from variation.schemas.classification_response_schema import (
    Classification, ClassificationType, AmplificationClassification
)
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import Validator


class Amplification(Validator):
    """The Insertion Validator Base class."""

    async def get_valid_invalid_results(
        self, classification: AmplificationClassification, transcripts: List
    ) -> List[ValidationResult]:
        # Does not require any validation
        return [ValidationResult(
            accession=None,
            classification=classification,
            is_valid=True,
            errors=[]
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
        self, classification: Classification, errors: List
    ) -> List:
        """Return empty list since amplification does not require transcripts

        :param Classification classification: A classification for a list of tokens
        :param List errors: List of errors
        :return: Empty list
        """
        return []
