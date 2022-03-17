"""Module for Protein Termination classification."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from .set_based_classifier import SetBasedClassifier


class ProteinTerminationClassifier(SetBasedClassifier):
    """The protein termination classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return Protein Termination classification type."""
        return ClassificationType.PROTEIN_TERMINATION

    def exact_match_candidates(self) -> List[List[str]]:
        """Return list of tokens for protein termination classification."""
        return [
            ["GeneSymbol", "ProteinTermination"],
            ["ProteinTermination"]
        ]
