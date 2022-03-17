"""Module for Genomic Substitution Translation."""
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import GenomicSubstitutionToken, Token


class GenomicSubstitution(Translator):
    """The Genomic Substitution Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Genomic Substitution."""
        return type == ClassificationType.GENOMIC_SUBSTITUTION

    def is_token_instance(self, token: Token) -> bool:
        """Return if the token is an Genomic Substitution token instance."""
        return isinstance(token, GenomicSubstitutionToken)
