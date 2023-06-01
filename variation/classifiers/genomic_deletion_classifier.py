"""A module for the Genomic Deletion Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import TokenType
from variation.classifiers import SetBasedClassifier


class GenomicDeletionClassifier(SetBasedClassifier):
    """The Genomic Deletion Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Deletion classification type."""
        return ClassificationType.GENOMIC_DELETION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.CHROMOSOME, TokenType.GENOMIC_DELETION],
            [TokenType.GENE, TokenType.PROTEIN_SUBSTITUTION, TokenType.GENOMIC_DELETION],  # noqa: E501
            [TokenType.GENOMIC_DELETION, TokenType.GENE],
            [TokenType.GENE, TokenType.GENOMIC_DELETION],
            [TokenType.HGVS, TokenType.GENOMIC_DELETION],
            [TokenType.REFERENCE_SEQUENCE, TokenType.GENOMIC_DELETION],
            [TokenType.LOCUS_REFERENCE_GENOMIC, TokenType.GENOMIC_DELETION]
        ]
