"""Module for Genomic Deletion Range Translation."""
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import \
    GenomicDeletionRangeToken, Token


class GenomicDeletionRange(Translator):
    """The Genomic Insertion Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Genomic Insertion."""
        return type == ClassificationType.GENOMIC_DELETION_RANGE

    def is_token_instance(self, token: Token) -> bool:
        """Return if the token is an Genomic Deletion Range token instance."""
        return isinstance(token, GenomicDeletionRangeToken)
