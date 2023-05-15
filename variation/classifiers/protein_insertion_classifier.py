"""A module for the Protein Insertion Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import TokenType
from variation.classifiers import SetBasedClassifier


class ProteinInsertionClassifier(SetBasedClassifier):
    """The Protein Insertion Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Protein Insertion classification type."""
        return ClassificationType.PROTEIN_INSERTION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.PROTEIN_INSERTION],
            [TokenType.GENE, TokenType.PROTEIN_SUBSTITUTION, TokenType.PROTEIN_INSERTION],  # noqa: E501
            [TokenType.PROTEIN_INSERTION, TokenType.GENE],
            [TokenType.GENE, TokenType.PROTEIN_INSERTION],
            [TokenType.HGVS, TokenType.PROTEIN_INSERTION],
            [TokenType.REFERENCE_SEQUENCE, TokenType.PROTEIN_INSERTION],
            [TokenType.LOCUS_REFERENCE_GENOMIC, TokenType.PROTEIN_INSERTION]
        ]
