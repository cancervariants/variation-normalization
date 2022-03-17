"""Module for Protein Substitution Translation."""
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import ProteinSubstitutionToken, Token


class ProteinSubstitution(Translator):
    """The Protein Substitution Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Protein Substitution."""
        return type == ClassificationType.PROTEIN_SUBSTITUTION

    def is_token_instance(self, token: Token) -> bool:
        """Return if the token is an Protein Substitution token instance."""
        return isinstance(token, ProteinSubstitutionToken)
