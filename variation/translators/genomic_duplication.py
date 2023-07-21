"""Module for Genomic Duplication Translation."""
from variation.translators.genomic_del_dup_base import GenomicDelDupTranslator
from variation.schemas.classification_response_schema import ClassificationType


class GenomicDuplication(GenomicDelDupTranslator):
    """The Genomic Duplication Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Genomic Duplication."""
        return type == ClassificationType.GENOMIC_DUPLICATION
