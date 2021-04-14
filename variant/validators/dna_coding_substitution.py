"""The module for Coding DNA Substitution Validation."""
from .dna_sequence_variant_base import DNASequenceVariantBase
from variant.schemas.classification_response_schema import \
    ClassificationType
from variant.schemas.token_response_schema import CodingDNASubstitutionToken


class CodingDNASubstitution(DNASequenceVariantBase):
    """The Coding DNA Substitution Validator class."""

    def variant_name(self):
        """Return the variant name."""
        return 'coding dna substitution'

    def is_token_instance(self, t):
        """Check that token is Coding DNA Substitution."""
        return t.token_type == 'CodingDNASubstitution'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is amino acid
        substitution.
        """
        return classification_type == ClassificationType.DNA_CODING_SUBSTITUTION  # noqa: E501

    def human_description(self, transcript,
                          psub_token: CodingDNASubstitutionToken) -> str:
        """Return a human description of the identified variant."""
        return f'An coding DNA substitution from {psub_token.ref_nucleotide}' \
               f' to {psub_token.new_nucleotide} at position ' \
               f'{psub_token.position} on transcript {transcript}'
