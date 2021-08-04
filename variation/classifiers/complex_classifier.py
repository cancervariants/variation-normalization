"""Module for Complex Classifier."""
from typing import List
from .set_based_classifier import SetBasedClassifier
from variation.schemas.classification_response_schema import ClassificationType


class ComplexClassifier(SetBasedClassifier):
    """The Complex Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the complex classification type."""
        return ClassificationType.COMPLEX

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the token matches for complex classification."""
        return [
          ['GeneSymbol', 'Amplification'],
          ['GeneSymbol', 'Exon'],
          ['GeneSymbol', 'Exon', 'Deletion'],
          ['Amplification'],
        ]
