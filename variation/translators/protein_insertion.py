"""Module for Protein Insertion Translation."""
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import ProteinInsertionToken, Token


class ProteinInsertion(Translator):
    """The Protein Insertion Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Protein Insertion."""
        return type == ClassificationType.PROTEIN_INSERTION

    def is_token_instance(self, token: Token) -> bool:
        """Return if the token is an Protein Insertion token instance."""
        return isinstance(token, ProteinInsertionToken)
