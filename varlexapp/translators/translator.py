from abc import ABC, abstractmethod

from ..models import ValidationResult, ClassificationType, VariantRepresentation

class Translator(ABC):

    @abstractmethod
    def translate(self, res: ValidationResult) -> VariantRepresentation:
        pass

    @abstractmethod
    def can_translate(self, type: ClassificationType) -> bool:
        pass
