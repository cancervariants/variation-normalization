"""Module for Coding DNA Reference Agree Translation."""
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import CodingDNAReferenceAgreeToken, Token


class CodingDNAReferenceAgree(Translator):
    """The Coding DNA Reference Agree Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Coding DNA Reference Agree."""
        return type == ClassificationType.CODING_DNA_REFERENCE_AGREE

    def is_token_instance(self, token: Token) -> bool:
        """Return if the token is an Coding DNA Reference Agree token
        instance.
        """
        return isinstance(token, CodingDNAReferenceAgreeToken)
