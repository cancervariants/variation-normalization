"""A module for the Amino Acid Deletion Classifier."""
from typing import List
from .set_based_classifier import SetBasedClassifier
from variation.schemas.classification_response_schema import ClassificationType


class AminoAcidDeletionClassifier(SetBasedClassifier):
    """The Amino Acid Deletion Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Amino Acid Deletion classification type."""
        return ClassificationType.AMINO_ACID_DELETION

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ['AminoAcidDeletion'],
            ['GeneSymbol', 'AminoAcidSubstitution', 'AminoAcidDeletion'],
            ['AminoAcidDeletion', 'GeneSymbol'],
            ['GeneSymbol', 'AminoAcidDeletion'],
            ['HGVS', 'AminoAcidDeletion'],
            ['ReferenceSequence', 'AminoAcidDeletion'],
            ['LocusReferenceGenomic', 'AminoAcidDeletion']
        ]
