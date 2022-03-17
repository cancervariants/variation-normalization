"""A module for the Silent Mutation Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from .set_based_classifier import SetBasedClassifier


class SilentMutationClassifier(SetBasedClassifier):
    """The Silent Mutation Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Silent Mutation classification type."""
        return ClassificationType.SILENT_MUTATION

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ["SilentMutation"],
            ["GeneSymbol", "SilentMutation"],
            ["HGVS", "SilentMutation"],
            ["ReferenceSequence", "SilentMutation"]
        ]
