"""Module for Validation."""
from typing import List
from abc import ABC, abstractmethod
from varlexapp.schemas.classification_response_schema import Classification, \
    ClassificationType
from varlexapp.schemas.validation_response_schema import ValidationResult


class Validator(ABC):
    """The validator class."""

    @abstractmethod
    def validate(self, classification: Classification) \
            -> List[ValidationResult]:
        """Return validation result for a given classification."""
        pass

    @abstractmethod
    def validates_classification_type(self,
                                      classification_type: ClassificationType)\
            -> bool:
        """Check that classification type matches."""
        pass
