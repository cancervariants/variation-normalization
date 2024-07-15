"""Module for Classification methods."""

from abc import ABC, abstractmethod

from variation.schemas.classification_response_schema import Classification
from variation.schemas.token_response_schema import Token, TokenType


class Classifier(ABC):
    """The Classifier class."""

    @abstractmethod
    def match(self, tokens: list[Token]) -> Classification | None:
        """Return the classification from a list of token matches.

        :param tokens: List of ordered tokens that are exact match candidates for a
            given classification
        :return: A classification for the list of matched tokens
        """

    @abstractmethod
    def exact_match_candidates(self) -> list[list[TokenType]]:
        """Return the token match candidates for a given classification.

        :return: List of list of tokens, where order matters, that represent a given
            classification.
        """

    def can_classify(self, tokens: list[Token]) -> bool:
        """Return whether or not a list of tokens can be classified by a given
        classification

        :param tokens: List of tokens found in an input query
        :return: `True` if a list of tokens matches the tokens needed, where order
            matters, to represent a given classification. `False`, otherwise.
        """
        token_types = [t.token_type for t in tokens]
        exact_matches: list[list[TokenType]] = [
            c for c in self.exact_match_candidates() if token_types == c
        ]

        return len(exact_matches) == 1
