"""Module for Amino Acid Substitution Translation."""
from variant.translators.translator import Translator
from variant.schemas.classification_response_schema import ClassificationType
from variant.schemas.token_response_schema import AminoAcidSubstitutionToken


class AminoAcidSubstitution(Translator):
    """The Amino Acid Substitution Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Amino Acid Substitution."""
        return type == ClassificationType.AMINO_ACID_SUBSTITUTION

    def is_token_instance(self, token):
        """Return if the token is an Amino Acid Substitution token instance."""
        return isinstance(token, AminoAcidSubstitutionToken)
