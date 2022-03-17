"""A module for the Protein Deletion Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from .set_based_classifier import SetBasedClassifier


class ProteinDeletionClassifier(SetBasedClassifier):
    """The Protein Deletion Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Protein Deletion classification type."""
        return ClassificationType.PROTEIN_DELETION

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ["ProteinDeletion"],
            ["GeneSymbol", "ProteinSubstitution", "ProteinDeletion"],
            ["ProteinDeletion", "GeneSymbol"],
            ["GeneSymbol", "ProteinDeletion"],
            ["HGVS", "ProteinDeletion"],
            ["ReferenceSequence", "ProteinDeletion"],
            ["LocusReferenceGenomic", "ProteinDeletion"]
        ]
