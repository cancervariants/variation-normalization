"""Module for Genomic Deletion Translation."""
from variant.translators.translator import Translator
from variant.schemas.classification_response_schema import ClassificationType
from variant.schemas.token_response_schema import GenomicDeletionToken


class GenomicDeletion(Translator):
    """The Genomic Deletion Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Genomic Deletion."""
        return type == ClassificationType.GENOMIC_DELETION

    def is_token_instance(self, token):
        """Return if the token is an Genomic Deletion token instance."""
        return isinstance(token, GenomicDeletionToken)
