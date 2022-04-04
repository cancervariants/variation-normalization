"""Module for Genomic Uncertain Deletion Translation."""
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import\
    GenomicUncertainDeletionToken, Token


class GenomicUncertainDeletion(Translator):
    """The Genomic Uncertain Deletion Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Genomic Uncertain Deletion."""
        return type == ClassificationType.GENOMIC_UNCERTAIN_DELETION

    def is_token_instance(self, token: Token) -> bool:
        """Return if the token is an Genomic Uncertain Deletion token
        instance.
        """
        return isinstance(token, GenomicUncertainDeletionToken)
