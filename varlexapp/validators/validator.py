from typing import List
from abc import ABC, abstractmethod

from ..models import ValidationResult, Classification, ClassificationType

class Validator(ABC):
    @abstractmethod
    def validate(self, classification: Classification) -> List[ValidationResult]:
        pass

    @abstractmethod
    def validates_classification_type(self, classification_type: ClassificationType) -> bool:
        pass
