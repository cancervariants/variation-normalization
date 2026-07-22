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

    @staticmethod
    def split_gene_cdna_tokens(tokens: list[Token]) -> tuple[Token, Token]:
        """Split parsed HGVS tokens into gene and cDNA change token components.

        Supports token lists with either two elements
        (``gene``, ``cdna_change``) or three elements
        (``gene``, ``protein_change``, ``cdna_change``).

        :param tokens: Parsed HGVS tokens
        :return: A tuple containing the gene symbol and cDNA change tokens
        """
        if len(tokens) == 2:
            gene_token, cdna_change_token = tokens
        else:
            gene_token, _, cdna_change_token = tokens

        return gene_token, cdna_change_token

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
