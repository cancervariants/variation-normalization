"""A module for the Coding DNA Substitution Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from .set_based_classifier import SetBasedClassifier


class CodingDNASubstitutionClassifier(SetBasedClassifier):
    """The Coding DNA Substitution Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Coding DNA Substitution classification type."""
        return ClassificationType.CODING_DNA_SUBSTITUTION

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ["CodingDNASubstitution"],
            ["GeneSymbol", "ProteinSubstitution", "CodingDNASubstitution"],
            ["CodingDNASubstitution", "GeneSymbol"],
            ["GeneSymbol", "CodingDNASubstitution"],
            ["HGVS", "CodingDNASubstitution"],
            ["ReferenceSequence", "CodingDNASubstitution"],
            ["LocusReferenceGenomic", "CodingDNASubstitution"]
        ]
