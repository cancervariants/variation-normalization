"""A module for the Genomic Uncertain Deletion Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from .set_based_classifier import SetBasedClassifier


class GenomicUncertainDeletionClassifier(SetBasedClassifier):
    """The Genomic Uncertain Deletion Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Uncertain Deletion classification type."""
        return ClassificationType.GENOMIC_UNCERTAIN_DELETION

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ["Chromosome", "GenomicUncertainDeletion"],
            ["GeneSymbol", "GenomicUncertainDeletion"],
            ["GenomicUncertainDeletion", "GeneSymbol"],
            ["HGVS", "GenomicUncertainDeletion"],
            ["ReferenceSequence", "GenomicUncertainDeletion"],
            ["LocusReferenceGenomic", "GenomicUncertainDeletion"]
        ]
