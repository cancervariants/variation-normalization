"""Module for Silent Mutation Translation."""
from .polypeptide_sequence_variant_base import PolypeptideSequenceVariantBase
from variant.schemas.classification_response_schema import ClassificationType
from variant.schemas.token_response_schema import SilentMutationToken


class SilentMutation(PolypeptideSequenceVariantBase):
    """The Amino Acid Substitution Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Silent Mutation."""
        return type == ClassificationType.SILENT_MUTATION

    def is_token_instance(self, token):
        """Return if the token is an Silent Mutation token instance."""
        return isinstance(token, SilentMutationToken)
