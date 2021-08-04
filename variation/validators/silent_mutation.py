"""The module for Silent Mutation Validation."""
from .polypeptide_sequence_variation_base import PolypeptideSequenceVariationBase
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import SilentMutationToken


class SilentMutation(PolypeptideSequenceVariationBase):
    """The Silent Mutation Validator class."""

    def variation_name(self):
        """Return the variation name."""
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
        """Return a human description of the identified variation."""
        return f'A silent mutation from {psub_token.ref_protein}' \
               f' to {psub_token.alt_protein} at position ' \
               f'{psub_token.position} on transcript {transcript}'
