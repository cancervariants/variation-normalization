"""The module for Silent Mutation Validation."""
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import Token
from .polypeptide_sequence_variation_base import PolypeptideSequenceVariationBase


class SilentMutation(PolypeptideSequenceVariationBase):
    """The Silent Mutation Validator class."""

    def variation_name(self) -> str:
        """Return the variation name."""
        return "silent mutation"

    def is_token_instance(self, t: Token) -> bool:
        """Check that token is Silent Mutation."""
        return t.token_type == "SilentMutation"

    def validates_classification_type(
            self, classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is silent
        mutation.
        """
        return classification_type == ClassificationType.SILENT_MUTATION
