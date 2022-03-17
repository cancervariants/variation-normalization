"""Module for classifying Protein DelIns."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from .set_based_classifier import SetBasedClassifier


class ProteinDelinsClassifier(SetBasedClassifier):
    """The Protein DeslIns Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Protein DelIns classification type."""
        return ClassificationType.PROTEIN_DELINS

    def exact_match_candidates(self) -> List[List[str]]:
        """Return Protein DelIns token combinations."""
        return [
            ["GeneSymbol", "ProteinDelins"],
            ["ProteinDelins"]
        ]
