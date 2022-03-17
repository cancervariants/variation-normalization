"""Module for Polypeptide Truncation Translation."""
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import PolypeptideTruncationToken, Token


class PolypeptideTruncation(Translator):
    """The Polypeptide Truncation Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Polypeptide Truncation."""
        return type == ClassificationType.POLYPEPTIDE_TRUNCATION

    def is_token_instance(self, token: Token) -> bool:
        """Return if the token is an Polypeptide Truncation token instance."""
        return isinstance(token, PolypeptideTruncationToken)
