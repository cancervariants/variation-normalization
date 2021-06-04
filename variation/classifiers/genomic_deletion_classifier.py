"""A module for the Genomic Deletion Classifier."""
from typing import List
from .set_based_classifier import SetBasedClassifier
from variation.schemas.classification_response_schema import ClassificationType


class GenomicDeletionClassifier(SetBasedClassifier):
    """The Genomic Deletion Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Deletion classification type."""
        return ClassificationType.GENOMIC_DELETION

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ['GenomicDeletion'],
            ['GeneSymbol', 'AminoAcidSubstitution', 'GenomicDeletion'],
            ['GenomicDeletion', 'GeneSymbol'],
            ['GeneSymbol', 'GenomicDeletion'],
            ['HGVS', 'GenomicDeletion'],
            ['ReferenceSequence', 'GenomicDeletion'],
            ['LocusReferenceGenomic', 'GenomicDeletion']
        ]
