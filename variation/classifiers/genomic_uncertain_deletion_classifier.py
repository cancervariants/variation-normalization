"""A module for the Genomic Uncertain Deletion Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import TokenType
from variation.classifiers import Classifier


class GenomicUncertainDeletionClassifier(Classifier):
    """The Genomic Uncertain Deletion Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Uncertain Deletion classification type."""
        return ClassificationType.GENOMIC_UNCERTAIN_DELETION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.GENOMIC_UNCERTAIN_DELETION],
            [TokenType.GENOMIC_UNCERTAIN_DELETION, TokenType.GENE],
            [TokenType.HGVS, TokenType.GENOMIC_UNCERTAIN_DELETION]
        ]
