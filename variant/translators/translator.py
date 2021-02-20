"""Module for translation."""
from abc import ABC, abstractmethod
from variant.schemas.ga4gh_vrs import Allele
from variant.schemas.validation_response_schema import ValidationResult
from variant.schemas.classification_response_schema import ClassificationType


class Translator(ABC):
    """The translation class."""

    @abstractmethod
    def translate(self, res: ValidationResult) -> Allele:
        """Translate a validation result to a VRS representation."""
        raise NotImplementedError

    @abstractmethod
    def can_translate(self, type: ClassificationType) -> bool:
        """Determine if it's possible to translate a classification."""
        raise NotImplementedError
