"""Module for translation."""
from abc import ABC, abstractmethod
from ..models import ValidationResult, VariantRepresentation
from varlexapp.schemas.classification_response_schema import ClassificationType


class Translator(ABC):
    """The translation class."""

    @abstractmethod
    def translate(self, res: ValidationResult) -> VariantRepresentation:
        """Translate a validation result to a VRS representation."""
        pass

    @abstractmethod
    def can_translate(self, type: ClassificationType) -> bool:
        """Determine if it's possible to translate a classification."""
        pass
