"""A module for the Protein DelIns Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from .set_based_classifier import SetBasedClassifier


class ProteinDelInsClassifier(SetBasedClassifier):
    """The Protein DelIns Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Protein DelIns classification type."""
        return ClassificationType.PROTEIN_DELINS

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ["ProteinDelIns"],
            ["GeneSymbol", "ProteinSubstitution", "ProteinDelIns"],
            ["ProteinDelIns", "GeneSymbol"],
            ["GeneSymbol", "ProteinDelIns"],
            ["HGVS", "ProteinDelIns"],
            ["ReferenceSequence", "ProteinDelIns"],
            ["LocusReferenceGenomic", "ProteinDelIns"]
        ]
