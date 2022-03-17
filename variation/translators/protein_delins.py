"""Module for Protein DelIns Translation."""
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import ProteinDelInsToken, Token


class ProteinDelIns(Translator):
    """The Protein DelIns Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Protein DelIns."""
        return type == ClassificationType.PROTEIN_DELINS

    def is_token_instance(self, token: Token) -> bool:
        """Return if the token is an Protein DelIns token instance."""
        return isinstance(token, ProteinDelInsToken)
