"""Module for Coding DNA DelIns Translation."""
from variant.translators.translator import Translator
from variant.schemas.classification_response_schema import ClassificationType
from variant.schemas.token_response_schema import CodingDNADelInsToken


class CodingDNADelins(Translator):
    """The Coding DNA Substitution Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Coding DNA DelIns."""
        return type == ClassificationType.CODING_DNA_DELINS

    def is_token_instance(self, token):
        """Return if the token is an Coding DNA DelIns token instance."""
        return isinstance(token, CodingDNADelInsToken)
