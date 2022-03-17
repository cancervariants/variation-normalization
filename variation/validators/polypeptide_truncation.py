"""The module for Polypeptide Truncation Validation."""
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import Token
from .polypeptide_sequence_variation_base import PolypeptideSequenceVariationBase


class PolypeptideTruncation(PolypeptideSequenceVariationBase):
    """The Polypeptide Truncation Validator class."""

    def variation_name(self) -> str:
        """Return the variation name."""
        return "polypeptide truncation"

    def is_token_instance(self, t: Token) -> bool:
        """Check that token is Polypeptide Truncation"""
        return t.token_type == "PolypeptideTruncation"

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is polypeptide
        truncation.
        """
        return classification_type == ClassificationType.POLYPEPTIDE_TRUNCATION
