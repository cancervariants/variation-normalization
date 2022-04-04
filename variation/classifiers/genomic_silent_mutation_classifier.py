"""A module for the Genomic Silent Mutation Classifier."""
from typing import List

from variation.schemas.classification_response_schema import ClassificationType
from .set_based_classifier import SetBasedClassifier


class GenomicSilentMutationClassifier(SetBasedClassifier):
    """The Genomic Silent Mutation Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Genomic Silent Mutation classification type."""
        return ClassificationType.GENOMIC_SILENT_MUTATION

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ["Chromosome", "GenomicSilentMutation"],
            ["GeneSymbol", "ProteinSubstitution", "GenomicSilentMutation"],
            ["GenomicSilentMutation", "GeneSymbol"],
            ["GeneSymbol", "GenomicSilentMutation"],
            ["HGVS", "GenomicSilentMutation"],
            ["ReferenceSequence", "GenomicSilentMutation"]
        ]
