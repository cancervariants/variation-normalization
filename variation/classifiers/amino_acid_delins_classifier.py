"""A module for the AminoAcid DelIns Classifier."""
from typing import List
from .set_based_classifier import SetBasedClassifier
from variation.schemas.classification_response_schema import ClassificationType


class AminoAcidDelInsClassifier(SetBasedClassifier):
    """The Amino Acid DelIns Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Amino Acid DelIns classification type."""
        return ClassificationType.AMINO_ACID_DELINS

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ['AminoAcidDelIns'],
            ['GeneSymbol', 'AminoAcidSubstitution', 'AminoAcidDelIns'],
            ['AminoAcidDelIns', 'GeneSymbol'],
            ['GeneSymbol', 'AminoAcidDelIns'],
            ['HGVS', 'AminoAcidDelIns'],
            ['ReferenceSequence', 'AminoAcidDelIns'],
            ['LocusReferenceGenomic', 'AminoAcidDelIns']
        ]
