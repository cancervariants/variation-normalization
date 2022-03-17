"""A module for the Genomic Insertion Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from .set_based_classifier import SetBasedClassifier


class GenomicInsertionClassifier(SetBasedClassifier):
    """The Genomic Insertion Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Insertion classification type."""
        return ClassificationType.GENOMIC_INSERTION

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ["Chromosome", "GenomicInsertion"],
            ["GeneSymbol", "ProteinSubstitution", "GenomicInsertion"],
            ["GenomicInsertion", "GeneSymbol"],
            ["GeneSymbol", "GenomicInsertion"],
            ["HGVS", "GenomicInsertion"],
            ["ReferenceSequence", "GenomicInsertion"],
            ["LocusReferenceGenomic", "GenomicInsertion"]
        ]
