"""A module for the Genomic Duplication Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from .set_based_classifier import SetBasedClassifier


class GenomicDuplicationClassifier(SetBasedClassifier):
    """The Genomic Duplication Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Duplication classification type."""
        return ClassificationType.GENOMIC_DUPLICATION

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ["Chromosome", "GenomicDuplication"],
            ["GenomicDuplication", "GeneSymbol"],
            ["GeneSymbol", "GenomicDuplication"],
            ["HGVS", "GenomicDuplication"],
            ["ReferenceSequence", "GenomicDuplication"],
            ["LocusReferenceGenomic", "GenomicDuplication"],
            ["GenomicDuplicationRange", "GeneSymbol"],
            ["GeneSymbol", "GenomicDuplicationRange"],
            ["HGVS", "GenomicDuplicationRange"],
            ["ReferenceSequence", "GenomicDuplicationRange"],
            ["LocusReferenceGenomic", "GenomicDuplicationRange"]
        ]
