"""A module for the Genomic Copy Number Loss Classifier."""
from typing import List
from .set_based_classifier import SetBasedClassifier
from variation.schemas.classification_response_schema import ClassificationType


class GenomicCopyNumberLossClassifier(SetBasedClassifier):
    """The Genomic Copy Number Loss Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Copy Number Loss classification type."""
        return ClassificationType.GENOMIC_COPY_NUMBER_LOSS

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ['GeneSymbol', 'GenomicCopyNumberLoss'],
            ['GenomicCopyNumberLoss', 'GeneSymbol'],
            ['HGVS', 'GenomicCopyNumberLoss'],
            ['ReferenceSequence', 'GenomicCopyNumberLoss'],
            ['LocusReferenceGenomic', 'GenomicCopyNumberLoss']
        ]
