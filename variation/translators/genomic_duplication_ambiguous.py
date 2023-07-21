"""Module for Genomic Duplication Ambiguous Translation."""
from variation.translators.ambiguous_translator_base import AmbiguousTranslator
from variation.schemas.classification_response_schema import ClassificationType


class GenomicDuplicationAmbiguous(AmbiguousTranslator):
    """The Genomic Duplication Ambiguous Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Genomic Duplication Ambiguous."""
        return type == ClassificationType.GENOMIC_DUPLICATION_AMBIGUOUS
