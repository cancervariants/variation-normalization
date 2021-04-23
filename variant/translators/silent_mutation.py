"""Module for Silent Mutation Translation."""
from variant.translators.translator import Translator
from variant.schemas.classification_response_schema import ClassificationType
from variant.schemas.token_response_schema import SilentMutationToken


class SilentMutation(Translator):
    """The Silent Mutation Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Silent Mutation."""
        return type == ClassificationType.SILENT_MUTATION

    def is_token_instance(self, token):
        """Return if the token is an Silent Mutation token instance."""
        return isinstance(token, SilentMutationToken)
