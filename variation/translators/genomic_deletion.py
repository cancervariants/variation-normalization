"""Module for Genomic Deletion Translation."""
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import GenomicDeletionToken, Token


class GenomicDeletion(Translator):
    """The Genomic Deletion Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Genomic Deletion."""
        return type == ClassificationType.GENOMIC_DELETION

    def is_token_instance(self, token: Token) -> bool:
        """Return if the token is an Genomic Deletion token instance."""
        return isinstance(token, GenomicDeletionToken)
