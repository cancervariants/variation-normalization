"""Module for Coding DNA Insertion Translation."""
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import CodingDNAInsertionToken


class CodingDNAInsertion(Translator):
    """The Coding DNA Insertion Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Coding DNA Insertion."""
        return type == ClassificationType.CODING_DNA_INSERTION

    def is_token_instance(self, token):
        """Return if the token is an Coding DNA Insertion token instance."""
        return isinstance(token, CodingDNAInsertionToken)
