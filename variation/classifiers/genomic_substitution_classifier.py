"""A module for the Genomic Substitution Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from .set_based_classifier import SetBasedClassifier


class GenomicSubstitutionClassifier(SetBasedClassifier):
    """The Genomic Substitution Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Substitution classification type."""
        return ClassificationType.GENOMIC_SUBSTITUTION

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ["Chromosome", "GenomicSubstitution"],
            ["GeneSymbol", "ProteinSubstitution", "GenomicSubstitution"],
            ["GenomicSubstitution", "GeneSymbol"],
            ["GeneSymbol", "GenomicSubstitution"],
            ["HGVS", "GenomicSubstitution"],
            ["ReferenceSequence", "GenomicSubstitution"]
        ]
