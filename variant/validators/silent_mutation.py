"""The module for Silent Mutation Validation."""
from .polypeptide_sequence_variant_base import PolypeptideSequenceVariantBase
from variant.schemas.classification_response_schema import \
    ClassificationType
from variant.schemas.token_response_schema import SilentMutationToken


class SilentMutation(PolypeptideSequenceVariantBase):
    """The Silent Mutation Validator class."""

    def variant_name(self):
        """Return the variant name."""
        return 'silent mutation'

    def is_token_instance(self, t):
        """Check that token is Silent Mutation."""
        return t.token_type == 'SilentMutation'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is silent
        mutation.
        """
        return classification_type == ClassificationType.SILENT_MUTATION

    def human_description(self, transcript,
                          psub_token: SilentMutationToken) -> str:
        """Return a human description of the identified variant."""
        return f'A silent mutation from {psub_token.ref_protein}' \
               f' to {psub_token.alt_protein} at position ' \
               f'{psub_token.position} on transcript {transcript}'
