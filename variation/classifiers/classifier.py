"""Module for Classification methods."""
from abc import ABC, abstractmethod
from typing import List, Optional

from variation.schemas.classification_response_schema import Classification
from variation.schemas.token_response_schema import Token, TokenType


class Classifier(ABC):
    """The Classifier class."""

    @abstractmethod
    def match(self, tokens: List[Token]) -> Optional[Classification]:
        """Return the classification from a list of token matches.

        :param tokens: List of ordered tokens that are exact match candidates for a
            given classification
        :return: A classification for the list of matched tokens
        """

    @abstractmethod
    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the token match candidates for a given classification.

        :return: List of list of tokens, where order matters, that represent a given
            classification.
        """
        pass

    def can_classify(self, tokens: List[Token]) -> bool:
        """Return whether or not a list of tokens can be classified by a given
        classification

        :param tokens: List of tokens found in an input query
        :return: `True` if a list of tokens matches the tokens needed, where order
            matters, to represent a given classification. `False`, otherwise.
        """
        token_types = list(map(lambda t: t.token_type, tokens))
        exact_matches: List[List[str]] = []

        for candidate in self.exact_match_candidates():
            if token_types == candidate:
                exact_matches.append(candidate)

        return len(exact_matches) == 1
