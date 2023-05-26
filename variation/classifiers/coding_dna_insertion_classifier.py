"""A module for the Coding DNA insertion Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import TokenType
from variation.classifiers import Classifier


class CodingDNAInsertionClassifier(Classifier):
    """The Coding DNA insertion Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Coding DNA insertion classification type."""
        return ClassificationType.CODING_DNA_INSERTION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.CODING_DNA_INSERTION],
            [TokenType.GENE, TokenType.PROTEIN_SUBSTITUTION, TokenType.CODING_DNA_INSERTION],  # noqa: E501
            [TokenType.CODING_DNA_INSERTION, TokenType.GENE],
            [TokenType.GENE, TokenType.CODING_DNA_INSERTION],
            [TokenType.HGVS, TokenType.CODING_DNA_INSERTION]
        ]
