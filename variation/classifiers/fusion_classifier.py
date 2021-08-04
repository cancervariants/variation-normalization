"""Module for the Fusion Classifier."""
from typing import List
from .set_based_classifier import SetBasedClassifier
from variation.schemas.classification_response_schema import ClassificationType


class FusionClassifier(SetBasedClassifier):
    """The Fusion classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Fusion classification type."""
        return ClassificationType.FUSION

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the tokens that classify as Fusion."""
        return [
          ['GenePair', 'Fusion'],
          ['GenePair'],
          ['GeneSymbol', 'Fusion']
        ]
