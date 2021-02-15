"""Module for Classifier."""
from abc import ABC, abstractmethod
from typing import List, Optional
from ..models import Classification, ClassificationType
from varlexapp.schemas.token_response_schema import Token


class Classifier(ABC):
    """The Classifier class."""

    @abstractmethod
    def match(self, tokens: List[Token]) -> Optional[Classification]:
        """Return the classification from a list of token matches."""
        pass

    @abstractmethod
    def classification_type(self) -> ClassificationType:
        """Return the classification type."""
        pass
