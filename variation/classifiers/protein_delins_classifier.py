"""A module for the Protein DelIns Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import TokenType
from variation.classifiers import SetBasedClassifier


class ProteinDelInsClassifier(SetBasedClassifier):
    """The Protein DelIns Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Protein DelIns classification type."""
        return ClassificationType.PROTEIN_DELINS

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.PROTEIN_DELINS],
            [TokenType.GENE, TokenType.PROTEIN_SUBSTITUTION, TokenType.PROTEIN_DELINS],
            [TokenType.PROTEIN_DELINS, TokenType.GENE],
            [TokenType.GENE, TokenType.PROTEIN_DELINS],
            [TokenType.HGVS, TokenType.PROTEIN_DELINS],
            [TokenType.REFERENCE_SEQUENCE, TokenType.PROTEIN_DELINS],
            [TokenType.LOCUS_REFERENCE_GENOMIC, TokenType.PROTEIN_DELINS]
        ]
