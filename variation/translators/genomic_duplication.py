"""Module for Genomic Duplication Translation."""
from variation.schemas.classification_response_schema import ClassificationType
from variation.translators.genomic_del_dup_base import GenomicDelDupTranslator


class GenomicDuplication(GenomicDelDupTranslator):
    """The Genomic Duplication Translator class."""

    def can_translate(self, classification_type: ClassificationType) -> bool:
        """Determine if it's possible to translate a classification.

        :param classification_type: Classification type found
        :return: `True` if `classification_type` matches translator's classification
            type. Otherwise, `False`
        """
        return classification_type == ClassificationType.GENOMIC_DUPLICATION
