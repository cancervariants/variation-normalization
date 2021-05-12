"""Module for Coding DNA Deletion Translation."""
from variant.translators.translator import Translator
from variant.schemas.classification_response_schema import ClassificationType
from variant.schemas.token_response_schema import CodingDNADeletionToken


class CodingDNADeletion(Translator):
    """The Coding DNA Deletion Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Coding DNA Deletion."""
        return type == ClassificationType.CODING_DNA_DELETION

    def is_token_instance(self, token):
        """Return if the token is an Coding DNA Deletion token instance."""
        return isinstance(token, CodingDNADeletionToken)
