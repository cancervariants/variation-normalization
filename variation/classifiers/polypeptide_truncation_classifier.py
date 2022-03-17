"""A module for the Polypeptide Truncation Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from .set_based_classifier import SetBasedClassifier


class PolypeptideTruncationClassifier(SetBasedClassifier):
    """The Polypeptide Truncation Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Polypeptide Truncation classification type."""
        return ClassificationType.POLYPEPTIDE_TRUNCATION

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ["PolypeptideTruncation"],
            ["GeneSymbol", "PolypeptideTruncation"],
            ["HGVS", "PolypeptideTruncation"],
            ["ReferenceSequence", "PolypeptideTruncation"]
        ]
