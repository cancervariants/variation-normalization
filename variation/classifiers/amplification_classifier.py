"""A module for the Amplification Classifier"""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, AmplificationClassification, Nomenclature
)
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier


class AmplificationClassifier(Classifier):
    """The Amplification Classifier class"""

    def classification_type(self) -> ClassificationType:
        """Return the Amplification classification type"""
        return ClassificationType.AMPLIFICATION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the token match candidates for the amplification classification.

        :return: List of list of tokens, where order matters, that represent an
        amplification classification.
        """
        return [
            [TokenType.GENE, TokenType.AMPLIFICATION]
        ]

    def match(self, tokens: List[Token]) -> AmplificationClassification:
        """Return the amplification classification from a list of token matches.

        :param tokens: List of ordered tokens that are exact match candidates for an
            amplification classification
        :return: amplification classification for the list of matched tokens
        """
        return AmplificationClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=tokens[0]
        )
