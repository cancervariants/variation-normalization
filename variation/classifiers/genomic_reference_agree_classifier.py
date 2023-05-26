"""A module for the Genomic Reference Agree Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import TokenType
from variation.classifiers import Classifier


class GenomicReferenceAgreeClassifier(Classifier):
    """The Genomic Reference Agree Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Reference Agree classification type."""
        return ClassificationType.GENOMIC_REFERENCE_AGREE

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GENE, TokenType.PROTEIN_SUBSTITUTION, TokenType.GENOMIC_REFERENCE_AGREE],  # noqa: E501
            [TokenType.GENOMIC_REFERENCE_AGREE, TokenType.GENE],
            [TokenType.GENE, TokenType.GENOMIC_REFERENCE_AGREE],
            [TokenType.HGVS, TokenType.GENOMIC_REFERENCE_AGREE]
        ]
