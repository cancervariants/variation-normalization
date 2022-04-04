"""Module for Coding DNA Deletion Translation."""
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import CodingDNADeletionToken, Token


class CodingDNADeletion(Translator):
    """The Coding DNA Deletion Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Coding DNA Deletion."""
        return type == ClassificationType.CODING_DNA_DELETION

    def is_token_instance(self, token: Token) -> bool:
        """Return if the token is an Coding DNA Deletion token instance."""
        return isinstance(token, CodingDNADeletionToken)
