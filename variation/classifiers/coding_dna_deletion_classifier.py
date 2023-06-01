"""A module for the Coding DNA Deletion Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import TokenType
from variation.classifiers import SetBasedClassifier


class CodingDNADeletionClassifier(SetBasedClassifier):
    """The Coding DNA Deletion Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Coding DNA Deletion classification type."""
        return ClassificationType.CODING_DNA_DELETION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.CODING_DNA_DELETION],
            [TokenType.GENE, TokenType.PROTEIN_SUBSTITUTION, TokenType.CODING_DNA_DELETION],  # noqa: E501
            [TokenType.CODING_DNA_DELETION, TokenType.GENE],
            [TokenType.GENE, TokenType.CODING_DNA_DELETION],
            [TokenType.HGVS, TokenType.CODING_DNA_DELETION],
            [TokenType.REFERENCE_SEQUENCE, TokenType.CODING_DNA_DELETION],
            [TokenType.LOCUS_REFERENCE_GENOMIC, TokenType.CODING_DNA_DELETION]
        ]
