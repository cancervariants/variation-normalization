"""A module for the DNA Coding Silent Mutation Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from .set_based_classifier import SetBasedClassifier


class CodingDNASilentMutationClassifier(SetBasedClassifier):
    """The Coding DNA Silent Mutation Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Coding DNA Silent Mutation classification type."""
        return ClassificationType.CODING_DNA_SILENT_MUTATION

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ["CodingDNASilentMutation"],
            ["GeneSymbol", "ProteinSubstitution", "CodingDNASilentMutation"],
            ["CodingDNASilentMutation", "GeneSymbol"],
            ["GeneSymbol", "CodingDNASilentMutation"],
            ["HGVS", "CodingDNASilentMutation"],
            ["ReferenceSequence", "CodingDNASilentMutation"],
            ["LocusReferenceGenomic", "CodingDNASilentMutation"]
        ]
