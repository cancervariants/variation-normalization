"""A module for the Genomic Deletion Range Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from .set_based_classifier import SetBasedClassifier


class GenomicDeletionRangeClassifier(SetBasedClassifier):
    """The Genomic Deletion Range Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Deletion Range classification type."""
        return ClassificationType.GENOMIC_DELETION_RANGE

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ["Chromosome", "GenomicDeletionRange"],
            ["GeneSymbol", "GenomicDeletionRange"],
            ["GenomicDeletionRange", "GeneSymbol"],
            ["HGVS", "GenomicDeletionRange"],
            ["ReferenceSequence", "GenomicDeletionRange"],
            ["LocusReferenceGenomic", "GenomicDeletionRange"]
        ]
