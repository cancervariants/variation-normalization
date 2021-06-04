"""The module for Amino Acid Substitution Validation."""
from .polypeptide_sequence_variation_base import PolypeptideSequenceVariationBase
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import AminoAcidSubstitutionToken


class AminoAcidSubstitution(PolypeptideSequenceVariationBase):
    """The Amino Acid Substitution Validator class."""

    def variation_name(self):
        """Return the variation name."""
        return 'amino acid substitution'

    def is_token_instance(self, t):
        """Check that token is Amino Acid Substitution."""
        return t.token_type == 'AminoAcidSubstitution'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is amino acid
        substitution.
        """
        return classification_type == ClassificationType.AMINO_ACID_SUBSTITUTION  # noqa: E501

    def human_description(self, transcript,
                          psub_token: AminoAcidSubstitutionToken) -> str:
        """Return a human description of the identified variation."""
        return f'An amino acid substitution from {psub_token.ref_protein}' \
               f' to {psub_token.alt_protein} at position ' \
               f'{psub_token.position} on transcript {transcript}'
