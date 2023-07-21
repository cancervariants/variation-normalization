"""Module for Genomic Deletion Ambiguous Translation."""
from variation.translators.ambiguous_translator_base import AmbiguousTranslator
from variation.schemas.classification_response_schema import ClassificationType


class GenomicDeletionAmbiguous(AmbiguousTranslator):
    """The Genomic Deletion Ambiguous Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Genomic Deletion Ambiguous."""
        return type == ClassificationType.GENOMIC_DELETION_AMBIGUOUS
