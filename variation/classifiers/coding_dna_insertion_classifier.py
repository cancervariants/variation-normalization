"""A module for the Coding DNA insertion Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from .set_based_classifier import SetBasedClassifier


class CodingDNAInsertionClassifier(SetBasedClassifier):
    """The Coding DNA insertion Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Coding DNA insertion classification type."""
        return ClassificationType.CODING_DNA_INSERTION

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ["CodingDNAInsertion"],
            ["GeneSymbol", "ProteinSubstitution", "CodingDNAInsertion"],
            ["CodingDNAInsertion", "GeneSymbol"],
            ["GeneSymbol", "CodingDNAInsertion"],
            ["HGVS", "CodingDNAInsertion"],
            ["ReferenceSequence", "CodingDNAInsertion"],
            ["LocusReferenceGenomic", "CodingDNAInsertion"]
        ]
