"""A module for the Amplification Classifier"""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from variation.classifiers import SetBasedClassifier


class AmplificationClassifier(SetBasedClassifier):
    """The Amplification Classifier class"""

    def classification_type(self) -> ClassificationType:
        """Return the Amplification classification type"""
        return ClassificationType.AMPLIFICATION

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ["GeneSymbol", "Amplification"]
        ]
