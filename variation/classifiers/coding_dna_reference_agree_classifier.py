"""A module for the DNA Coding Reference Agree Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import TokenType
from variation.classifiers import Classifier


class CodingDNAReferenceAgreeClassifier(Classifier):
    """The Coding DNA Reference Agree Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Coding DNA Reference Agree classification type."""
        return ClassificationType.CODING_DNA_REFERENCE_AGREE

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.CODING_DNA_REFERENCE_AGREE],
            [TokenType.GENE, TokenType.PROTEIN_SUBSTITUTION, TokenType.CODING_DNA_REFERENCE_AGREE],  # noqa: E501
            [TokenType.CODING_DNA_REFERENCE_AGREE, TokenType.GENE],
            [TokenType.GENE, TokenType.CODING_DNA_REFERENCE_AGREE],
            [TokenType.HGVS, TokenType.CODING_DNA_REFERENCE_AGREE]
        ]
