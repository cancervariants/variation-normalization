"""Module for Protein Alternate classification."""
from typing import List
from .set_based_classifier import SetBasedClassifier
from variation.schemas.classification_response_schema import ClassificationType


class ProteinAlternateClassifier(SetBasedClassifier):
    """The protein alternate classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Protein Alternate classification type."""
        return ClassificationType.PROTEIN_ALTERNATE

    def exact_match_candidates(self) -> List[List[str]]:
        """Return tokens that match Protein Alternate classification type."""
        return [
          ['GeneSymbol', 'ProteinAlternate'],
          ['ProteinAlternate']
        ]
