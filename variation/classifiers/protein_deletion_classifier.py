"""A module for the Protein Deletion Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import TokenType
from variation.classifiers import SetBasedClassifier


class ProteinDeletionClassifier(SetBasedClassifier):
    """The Protein Deletion Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Protein Deletion classification type."""
        return ClassificationType.PROTEIN_DELETION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.PROTEIN_DELETION],
            [TokenType.GENE, TokenType.PROTEIN_SUBSTITUTION, TokenType.PROTEIN_DELETION],  # noqa: E501
            [TokenType.PROTEIN_DELETION, TokenType.GENE],
            [TokenType.GENE, TokenType.PROTEIN_DELETION],
            [TokenType.HGVS, TokenType.PROTEIN_DELETION],
            [TokenType.REFERENCE_SEQUENCE, TokenType.PROTEIN_DELETION],
            [TokenType.LOCUS_REFERENCE_GENOMIC, TokenType.PROTEIN_DELETION]
        ]
