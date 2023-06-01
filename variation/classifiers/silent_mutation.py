"""A module for the Silent Mutation Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import TokenType
from variation.classifiers import SetBasedClassifier


class SilentMutationClassifier(SetBasedClassifier):
    """The Silent Mutation Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Silent Mutation classification type."""
        return ClassificationType.SILENT_MUTATION

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.SILENT_MUTATION],
            [TokenType.GENE, TokenType.SILENT_MUTATION],
            [TokenType.HGVS, TokenType.SILENT_MUTATION],
            [TokenType.REFERENCE_SEQUENCE, TokenType.SILENT_MUTATION]
        ]
