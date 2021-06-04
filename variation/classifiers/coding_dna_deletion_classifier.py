"""A module for the Coding DNA Deletion Classifier."""
from typing import List
from .set_based_classifier import SetBasedClassifier
from variation.schemas.classification_response_schema import ClassificationType


class CodingDNADeletionClassifier(SetBasedClassifier):
    """The Coding DNA Deletion Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Coding DNA Deletion classification type."""
        return ClassificationType.CODING_DNA_DELETION

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ['CodingDNADeletion'],
            ['GeneSymbol', 'AminoAcidSubstitution', 'CodingDNADeletion'],
            ['CodingDNADeletion', 'GeneSymbol'],
            ['GeneSymbol', 'CodingDNADeletion'],
            ['HGVS', 'CodingDNADeletion'],
            ['ReferenceSequence', 'CodingDNADeletion'],
            ['LocusReferenceGenomic', 'CodingDNADeletion']
        ]
