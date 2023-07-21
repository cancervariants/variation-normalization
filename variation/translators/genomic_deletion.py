"""Module for Genomic Deletion Translation."""
from variation.translators.genomic_del_dup_base import GenomicDelDupTranslator
from variation.schemas.classification_response_schema import ClassificationType


class GenomicDeletion(GenomicDelDupTranslator):
    """The Genomic Deletion Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Genomic Deletion."""
        return type == ClassificationType.GENOMIC_DELETION
