"""Module for Coding DNA Silent Mutation Translation."""
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import CodingDNASilentMutationToken, Token


class CodingDNASilentMutation(Translator):
    """The Coding DNA Silent Mutation Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Coding DNA Silent Mutation."""
        return type == ClassificationType.CODING_DNA_SILENT_MUTATION

    def is_token_instance(self, token: Token) -> bool:
        """Return if the token is an Coding DNA Silent Mutation token
        instance.
        """
        return isinstance(token, CodingDNASilentMutationToken)
