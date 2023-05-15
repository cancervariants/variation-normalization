"""A module for the Genomic Substitution Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import TokenType
from variation.classifiers import SetBasedClassifier


class GenomicSubstitutionClassifier(SetBasedClassifier):
    """The Genomic Substitution Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Substitution classification type."""
        return ClassificationType.GENOMIC_SUBSTITUTION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.CHROMOSOME, TokenType.GENOMIC_SUBSTITUTION],
            [TokenType.GENE, TokenType.PROTEIN_SUBSTITUTION, TokenType.GENOMIC_SUBSTITUTION],  # noqa: E501
            [TokenType.GENOMIC_SUBSTITUTION, TokenType.GENE],
            [TokenType.GENE, TokenType.GENOMIC_SUBSTITUTION],
            [TokenType.HGVS, TokenType.GENOMIC_SUBSTITUTION],
            [TokenType.REFERENCE_SEQUENCE, TokenType.GENOMIC_SUBSTITUTION]
        ]
