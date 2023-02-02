"""Module for Amplification Translation."""
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import AmplificationToken, Token


class Amplification(Translator):
    """The Amplification Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Amplification."""
        return type == ClassificationType.AMPLIFICATION

    def is_token_instance(self, token: Token) -> bool:
        """Return if the token is an Amplification token instance."""
        return isinstance(token, AmplificationToken)
