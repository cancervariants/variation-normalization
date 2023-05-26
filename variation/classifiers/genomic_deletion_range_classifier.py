"""A module for the Genomic Deletion Range Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import TokenType
from variation.classifiers import Classifier


class GenomicDeletionRangeClassifier(Classifier):
    """The Genomic Deletion Range Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Deletion Range classification type."""
        return ClassificationType.GENOMIC_DELETION_RANGE

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.GENOMIC_DELETION_RANGE],
            [TokenType.GENOMIC_DELETION_RANGE, TokenType.GENE],
            [TokenType.HGVS, TokenType.GENOMIC_DELETION_RANGE]
        ]
