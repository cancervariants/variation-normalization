"""A module for the Genomic Deletion Range Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import TokenType
from variation.classifiers import SetBasedClassifier


class GenomicDeletionRangeClassifier(SetBasedClassifier):
    """The Genomic Deletion Range Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Deletion Range classification type."""
        return ClassificationType.GENOMIC_DELETION_RANGE

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.CHROMOSOME, TokenType.GENOMIC_DELETION_RANGE],
            [TokenType.GENE, TokenType.GENOMIC_DELETION_RANGE],
            [TokenType.GENOMIC_DELETION_RANGE, TokenType.GENE],
            [TokenType.HGVS, TokenType.GENOMIC_DELETION_RANGE],
            [TokenType.REFERENCE_SEQUENCE, TokenType.GENOMIC_DELETION_RANGE],
            [TokenType.LOCUS_REFERENCE_GENOMIC, TokenType.GENOMIC_DELETION_RANGE]
        ]
