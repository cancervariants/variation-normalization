"""A module for the Amino Acid Insertion Classifier."""
from typing import List
from .set_based_classifier import SetBasedClassifier
from variation.schemas.classification_response_schema import ClassificationType


class AminoAcidInsertionClassifier(SetBasedClassifier):
    """The Amino Acid Insertion Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Amino Acid Insertion classification type."""
        return ClassificationType.AMINO_ACID_INSERTION

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ['AminoAcidInsertion'],
            ['GeneSymbol', 'AminoAcidSubstitution', 'AminoAcidInsertion'],
            ['AminoAcidInsertion', 'GeneSymbol'],
            ['GeneSymbol', 'AminoAcidInsertion'],
            ['HGVS', 'AminoAcidInsertion'],
            ['ReferenceSequence', 'AminoAcidInsertion'],
            ['LocusReferenceGenomic', 'AminoAcidInsertion']
        ]
