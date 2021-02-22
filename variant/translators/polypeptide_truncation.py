"""Module for Polypeptide Truncation Translation."""
from .polypeptide_sequence_variant_base import PolypeptideSequenceVariantBase
from variant.schemas.classification_response_schema import ClassificationType
from variant.schemas.token_response_schema import PolypeptideTruncationToken


class PolypeptideTruncation(PolypeptideSequenceVariantBase):
    """The Polypeptide Truncation Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Polypeptide Truncation."""
        return type == ClassificationType.POLYPEPTIDE_TRUNCATION

    def is_token_instance(self, token):
        """Return if the token is an Polypeptide Truncation token instance."""
        return isinstance(token, PolypeptideTruncationToken)
