"""Module for Genomic Duplication Translation."""
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import \
    GenomicDuplicationRangeToken, GenomicDuplicationToken, Token


class GenomicDuplication(Translator):
    """The Genomic Insertion Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Genomic Insertion."""
        return type == ClassificationType.GENOMIC_DUPLICATION

    def is_token_instance(self, token: Token) -> bool:
        """Return if the token is an Genomic Duplication token instance."""
        return isinstance(token, (GenomicDuplicationToken,
                                  GenomicDuplicationRangeToken))
