"""Module for Classifier."""
from abc import ABC, abstractmethod
from typing import List, Optional
from variant.schemas.classification_response_schema import Classification,\
    ClassificationType
from variant.schemas.token_response_schema import Token


class Classifier(ABC):
    """The Classifier class."""

    @abstractmethod
    def match(self, tokens: List[Token]) -> Optional[Classification]:
        """Return the classification from a list of token matches."""
        raise NotImplementedError

    @abstractmethod
    def classification_type(self) -> ClassificationType:
        """Return the classification type."""
        raise NotImplementedError
