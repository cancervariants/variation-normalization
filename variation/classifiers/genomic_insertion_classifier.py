"""A module for the Genomic Insertion Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import Classifier


class GenomicInsertionClassifier(Classifier):
    """The Genomic Insertion Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Insertion classification type."""
        return ClassificationType.GENOMIC_INSERTION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.PROTEIN_SUBSTITUTION, TokenType.GENOMIC_INSERTION],  # noqa: E501
            [TokenType.GENOMIC_INSERTION, TokenType.GENE],
            [TokenType.GENE, TokenType.GENOMIC_INSERTION],
            [TokenType.HGVS, TokenType.GENOMIC_INSERTION]
        ]
