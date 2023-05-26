"""A module for the Genomic Duplication Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import TokenType
from variation.classifiers import Classifier


class GenomicDuplicationClassifier(Classifier):
    """The Genomic Duplication Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Duplication classification type."""
        return ClassificationType.GENOMIC_DUPLICATION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENOMIC_DUPLICATION, TokenType.GENE],
            [TokenType.GENE, TokenType.GENOMIC_DUPLICATION],
            [TokenType.HGVS, TokenType.GENOMIC_DUPLICATION],
            [TokenType.GENOMIC_DUPLICATION_RANGE, TokenType.GENE],
            [TokenType.GENE, TokenType.GENOMIC_DUPLICATION_RANGE],
            [TokenType.HGVS, TokenType.GENOMIC_DUPLICATION_RANGE]
        ]
