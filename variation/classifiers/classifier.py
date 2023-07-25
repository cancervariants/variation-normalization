"""Module for Classification methods."""
from abc import ABC, abstractmethod
from typing import List, Optional

from variation.schemas.classification_response_schema import Classification
from variation.schemas.token_response_schema import Token


class Classifier(ABC):
    """The Classifier class."""

    @abstractmethod
    def match(self, tokens: List[Token]) -> Optional[Classification]:
        """Return the classification from a list of token matches."""

    @abstractmethod
    def exact_match_candidates(self) -> List[List[str]]:
        """Return the token match candidates for a given classification.

        :return: List of tokens, where order matters, that make up a given
            classification.
        """
        pass

    def can_classify(self, tokens: List[Token]) -> bool:
        token_types = list(map(lambda t: t.token_type, tokens))
        exact_matches: List[List[str]] = []

        for candidate in self.exact_match_candidates():
            if token_types == candidate:
                exact_matches.append(candidate)

        return len(exact_matches) == 1
