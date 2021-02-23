"""Module for Validation."""
from typing import List
from abc import ABC, abstractmethod
from variant.schemas.classification_response_schema import Classification, \
    ClassificationType
from variant.schemas.validation_response_schema import ValidationResult


class Validator(ABC):
    """The validator class."""

    @abstractmethod
    def validate(self, classification: Classification) \
            -> List[ValidationResult]:
        """Return validation result for a given classification."""
        raise NotImplementedError

    @abstractmethod
    def validates_classification_type(self,
                                      classification_type: ClassificationType)\
            -> bool:
        """Check that classification type matches."""
        raise NotImplementedError
