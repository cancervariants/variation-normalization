"""Module for Silent Mutation Translation."""
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import SilentMutationToken, Token


class SilentMutation(Translator):
    """The Silent Mutation Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Silent Mutation."""
        return type == ClassificationType.SILENT_MUTATION

    def is_token_instance(self, token: Token) -> bool:
        """Return if the token is an Silent Mutation token instance."""
        return isinstance(token, SilentMutationToken)
