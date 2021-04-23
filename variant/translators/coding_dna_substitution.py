"""Module for Coding DNA Substitution Translation."""
from variant.translators.translator import Translator
from variant.schemas.classification_response_schema import ClassificationType
from variant.schemas.token_response_schema import CodingDNASubstitutionToken


class CodingDNASubstitution(Translator):
    """The Coding DNA Substitution Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Coding DNA Substitution."""
        return type == ClassificationType.CODING_DNA_SUBSTITUTION

    def is_token_instance(self, token):
        """Return if the token is an Coding DNA Substitution token instance."""
        return isinstance(token, CodingDNASubstitutionToken)
