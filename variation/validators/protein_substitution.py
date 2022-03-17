"""The module for Protein Substitution Validation."""
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import Token
from .polypeptide_sequence_variation_base import PolypeptideSequenceVariationBase


class ProteinSubstitution(PolypeptideSequenceVariationBase):
    """The Protein Substitution Validator class."""

    def variation_name(self) -> str:
        """Return the variation name."""
        return "protein substitution"

    def is_token_instance(self, t: Token) -> bool:
        """Check that token is Protein Substitution."""
        return t.token_type == "ProteinSubstitution"

    def validates_classification_type(
            self, classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is protein
        substitution.
        """
        return classification_type == ClassificationType.PROTEIN_SUBSTITUTION
