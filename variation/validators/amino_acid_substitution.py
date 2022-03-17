"""The module for Amino Acid Substitution Validation."""
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import Token
from .polypeptide_sequence_variation_base import PolypeptideSequenceVariationBase


class AminoAcidSubstitution(PolypeptideSequenceVariationBase):
    """The Amino Acid Substitution Validator class."""

    def variation_name(self) -> str:
        """Return the variation name."""
        return "amino acid substitution"

    def is_token_instance(self, t: Token) -> bool:
        """Check that token is Amino Acid Substitution."""
        return t.token_type == "AminoAcidSubstitution"

    def validates_classification_type(
            self, classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is amino acid
        substitution.
        """
        return classification_type == ClassificationType.AMINO_ACID_SUBSTITUTION
