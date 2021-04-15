"""A module for the DNA Coding Substitution Classifier."""
from typing import List
from .set_based_classifier import SetBasedClassifier
from variant.schemas.classification_response_schema import ClassificationType


class DNACodingSubstitutionClassifier(SetBasedClassifier):
    """The DNA Coding Substitution Classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the DNA Coding Substitution classification type."""
        return ClassificationType.DNA_CODING_SUBSTITUTION

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the exact match token type candidates."""
        return [
            ['CodingDNASubstitution'],
            ['GeneSymbol', 'AminoAcidSubstitution', 'CodingDNASubstitution'],
            ['CodingDNASubstitution', 'GeneSymbol'],
            ['GeneSymbol', 'CodingDNASubstitution'],
            ['HGVS', 'CodingDNASubstitution'],
            ['ReferenceSequence', 'CodingDNASubstitution']
        ]
