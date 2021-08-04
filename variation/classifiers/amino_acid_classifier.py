"""A module for the Amino Acid Substitution Classifier."""
from typing import List
from .set_based_classifier import SetBasedClassifier
from variation.schemas.classification_response_schema import ClassificationType


class AminoAcidSubstitutionClassifier(SetBasedClassifier):
    """The Amino Acid Substitution Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the Amino Acid Substitution classification type."""
        return ClassificationType.AMINO_ACID_SUBSTITUTION

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ['AminoAcidSubstitution'],
            ['GeneSymbol', 'AminoAcidSubstitution'],
            ['HGVS', 'AminoAcidSubstitution'],
            ['ReferenceSequence', 'AminoAcidSubstitution']
        ]
