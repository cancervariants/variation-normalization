"""Module for Genomic Deletion Ambiguous Translation."""
from variation.schemas.classification_response_schema import ClassificationType
from variation.translators.ambiguous_translator_base import AmbiguousTranslator


class GenomicDeletionAmbiguous(AmbiguousTranslator):
    """The Genomic Deletion Ambiguous Translator class."""

    def can_translate(self, classification_type: ClassificationType) -> bool:
        """Determine if it's possible to translate a classification.

        :param classification_type: Classification type found
        :return: `True` if `classification_type` matches translator's classification
            type. Otherwise, `False`
        """
        return classification_type == ClassificationType.GENOMIC_DELETION_AMBIGUOUS
