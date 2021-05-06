"""Module for Amino Acid DelIns Translation."""
from variant.translators.translator import Translator
from variant.schemas.classification_response_schema import ClassificationType
from variant.schemas.token_response_schema import AminoAcidDelInsToken


class AminoAcidDelIns(Translator):
    """The Amino Acid DelIns Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Amino Acid DelIns."""
        return type == ClassificationType.AMINO_ACID_DELINS

    def is_token_instance(self, token):
        """Return if the token is an Amino Acid DelIns token instance."""
        return isinstance(token, AminoAcidDelInsToken)
