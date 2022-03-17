"""A module for the Protein Insertion Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from .set_based_classifier import SetBasedClassifier


class ProteinInsertionClassifier(SetBasedClassifier):
    """The Protein Insertion Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Protein Insertion classification type."""
        return ClassificationType.PROTEIN_INSERTION

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ["ProteinInsertion"],
            ["GeneSymbol", "ProteinSubstitution", "ProteinInsertion"],
            ["ProteinInsertion", "GeneSymbol"],
            ["GeneSymbol", "ProteinInsertion"],
            ["HGVS", "ProteinInsertion"],
            ["ReferenceSequence", "ProteinInsertion"],
            ["LocusReferenceGenomic", "ProteinInsertion"]
        ]
