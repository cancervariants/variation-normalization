"""A module for the Protein Substitution Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import TokenType
from variation.classifiers import SetBasedClassifier


class ProteinSubstitutionClassifier(SetBasedClassifier):
    """The ProteinSubstitution Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Protein Substitution classification type."""
        return ClassificationType.PROTEIN_SUBSTITUTION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.PROTEIN_SUBSTITUTION],
            [TokenType.GENE, TokenType.PROTEIN_SUBSTITUTION],
            [TokenType.HGVS, TokenType.PROTEIN_SUBSTITUTION],
            [TokenType.REFERENCE_SEQUENCE, TokenType.PROTEIN_SUBSTITUTION]
        ]
