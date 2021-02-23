"""The module for Polypeptide Truncation Validation."""
from .polypeptide_sequence_variant_base import PolypeptideSequenceVariantBase
from variant.schemas.classification_response_schema import \
    ClassificationType
from variant.schemas.token_response_schema import PolypeptideTruncationToken


class PolypeptideTruncation(PolypeptideSequenceVariantBase):
    """The Polypeptide Truncation Validator class."""

    def variant_name(self):
        """Return the variant name."""
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
        """Return a human description of the identified variant."""
        return f'A polypeptide truncation from {psub_token.ref_protein}' \
               f' to {psub_token.alt_protein} at position ' \
               f'{psub_token.position} on transcript {transcript}'
