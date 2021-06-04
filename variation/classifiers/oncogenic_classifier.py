"""Module for Oncogenic classification."""
from typing import List
from .set_based_classifier import SetBasedClassifier
from variation.schemas.classification_response_schema import ClassificationType


class OncogenicClassifier(SetBasedClassifier):
    """The oncogenic classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the oncogenic classification type."""
        return ClassificationType.ONCOGENIC

    def exact_match_candidates(self) -> List[List[str]]:
        """Return tokens that classify as ocongenic."""
        return [
          ['GeneSymbol', 'LossOfFunction'],
          ['GeneSymbol', 'GainOfFunction'],
          ['LossOfFunction'],
          ['GainOfFunction']
        ]
