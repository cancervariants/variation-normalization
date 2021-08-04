"""Module for classifying Protein Frameshift."""
from typing import List
from .set_based_classifier import SetBasedClassifier
from variation.schemas.classification_response_schema import ClassificationType


class ProteinFrameshiftClassifier(SetBasedClassifier):
    """The Protein Frameshift Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Protein Frameshift classification type."""
        return ClassificationType.PROTEIN_FRAMESHIFT

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the tokens that classify as Protein Frameshift."""
        return [
          ['GeneSymbol', 'ProteinFrameshift'],
          ['ProteinFrameshift']
        ]
