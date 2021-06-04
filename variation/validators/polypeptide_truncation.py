"""The module for Polypeptide Truncation Validation."""
from .polypeptide_sequence_variation_base import PolypeptideSequenceVariationBase
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import PolypeptideTruncationToken


class PolypeptideTruncation(PolypeptideSequenceVariationBase):
    """The Polypeptide Truncation Validator class."""

    def variation_name(self):
        """Return the variation name."""
        return 'polypeptide truncation'

    def is_token_instance(self, t):
        """Check that token is Polypeptide Truncation"""
        return t.token_type == 'PolypeptideTruncation'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is polypeptide
        truncation.
        """
        return classification_type == ClassificationType.POLYPEPTIDE_TRUNCATION

    def human_description(self, transcript,
                          psub_token: PolypeptideTruncationToken) -> str:
        """Return a human description of the identified variation."""
        return f'A polypeptide truncation from {psub_token.ref_protein}' \
               f' to {psub_token.alt_protein} at position ' \
               f'{psub_token.position} on transcript {transcript}'
