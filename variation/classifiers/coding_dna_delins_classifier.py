"""A module for the Coding DNA DelIns Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import TokenType
from variation.classifiers import SetBasedClassifier


class CodingDNADelInsClassifier(SetBasedClassifier):
    """The Coding DNA DelIns Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Coding DNA DelIns classification type."""
        return ClassificationType.CODING_DNA_DELINS

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.CODING_DNA_DELINS],
            [TokenType.GENE, TokenType.PROTEIN_SUBSTITUTION, TokenType.CODING_DNA_DELINS],  # noqa: E501
            [TokenType.CODING_DNA_DELINS, TokenType.GENE],
            [TokenType.GENE, TokenType.CODING_DNA_DELINS],
            [TokenType.HGVS, TokenType.CODING_DNA_DELINS],
            [TokenType.REFERENCE_SEQUENCE, TokenType.CODING_DNA_DELINS],
            [TokenType.LOCUS_REFERENCE_GENOMIC, TokenType.CODING_DNA_DELINS]
        ]
