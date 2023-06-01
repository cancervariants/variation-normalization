"""A module for the Genomic DelIns Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import TokenType
from variation.classifiers import SetBasedClassifier


class GenomicDelInsClassifier(SetBasedClassifier):
    """The Genomic DelIns Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic DelIns classification type."""
        return ClassificationType.GENOMIC_DELINS

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.CHROMOSOME, TokenType.GENOMIC_DELINS],
            [TokenType.GENE, TokenType.PROTEIN_SUBSTITUTION, TokenType.GENOMIC_DELINS],
            [TokenType.GENOMIC_DELINS, TokenType.GENE],
            [TokenType.GENE, TokenType.GENOMIC_DELINS],
            [TokenType.HGVS, TokenType.GENOMIC_DELINS],
            [TokenType.REFERENCE_SEQUENCE, TokenType.GENOMIC_DELINS]
        ]
