"""A module for the Protein Substitution Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from .set_based_classifier import SetBasedClassifier


class ProteinSubstitutionClassifier(SetBasedClassifier):
    """The ProteinSubstitution Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Protein Substitution classification type."""
        return ClassificationType.PROTEIN_SUBSTITUTION

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ["ProteinSubstitution"],
            ["GeneSymbol", "ProteinSubstitution"],
            ["HGVS", "ProteinSubstitution"],
            ["ReferenceSequence", "ProteinSubstitution"]
        ]
