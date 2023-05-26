"""Module for Genomic Reference Agree Translation."""
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import GenomicReferenceAgreeToken, Token


class GenomicReferenceAgree(Translator):
    """The Genomic Reference Agree Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Genomic Reference Agree."""
        return type == ClassificationType.GENOMIC_REFERENCE_AGREE

    def is_token_instance(self, token: Token) -> bool:
        """Return if the token is an Genomic Reference Agree token
        instance.
        """
        return isinstance(token, GenomicReferenceAgreeToken)
