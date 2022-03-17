"""Module for Protein Deletion Translation."""
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import ProteinDeletionToken, Token


class ProteinDeletion(Translator):
    """The Protein Deletion Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Protein Deletion."""
        return type == ClassificationType.PROTEIN_DELETION

    def is_token_instance(self, token: Token) -> bool:
        """Return if the token is an Protein Deletion token instance."""
        return isinstance(token, ProteinDeletionToken)
