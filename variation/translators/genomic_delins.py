"""Module for Genomic DelIns Translation."""
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import GenomicDelInsToken, Token


class GenomicDelIns(Translator):
    """The Genomic DelIns Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Genomic DelIns."""
        return type == ClassificationType.GENOMIC_DELINS

    def is_token_instance(self, token: Token) -> bool:
        """Return if the token is an Genomic DelIns token instance."""
        return isinstance(token, GenomicDelInsToken)
