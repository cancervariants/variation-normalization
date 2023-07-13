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
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.AMPLIFICATION]
        ]

    def match(self, tokens: List[Token]):
        return AmplificationClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=tokens[0]
        )
