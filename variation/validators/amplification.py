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
        self, classification: AmplificationClassification, accessions: List
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

    async def get_accessions(
        self, classification: Classification, errors: List
    ) -> List:
        """Return empty list since amplification does not require accessions

        :param classification: The classification for list of tokens
        :param errors: List of errors
        :return: Empty list
        """
        return []
