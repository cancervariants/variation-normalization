"""A module for the Coding DNA Substitution Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import TokenType
from variation.classifiers import SetBasedClassifier


class CodingDNASubstitutionClassifier(SetBasedClassifier):
    """The Coding DNA Substitution Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Coding DNA Substitution classification type."""
        return ClassificationType.CODING_DNA_SUBSTITUTION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.CODING_DNA_SUBSTITUTION],
            [TokenType.GENE, TokenType.PROTEIN_SUBSTITUTION, TokenType.CODING_DNA_SUBSTITUTION],  # noqa: E501
            [TokenType.CODING_DNA_SUBSTITUTION, TokenType.GENE],
            [TokenType.GENE, TokenType.CODING_DNA_SUBSTITUTION],
            [TokenType.HGVS, TokenType.CODING_DNA_SUBSTITUTION],
            [TokenType.REFERENCE_SEQUENCE, TokenType.CODING_DNA_SUBSTITUTION],
            [TokenType.LOCUS_REFERENCE_GENOMIC, TokenType.CODING_DNA_SUBSTITUTION]
        ]
