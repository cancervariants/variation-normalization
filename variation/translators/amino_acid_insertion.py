"""Module for Amino Acid Insertion Translation."""
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.token_response_schema import AminoAcidInsertionToken


class AminoAcidInsertion(Translator):
    """The Amino Acid Insertion Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Amino Acid Insertion."""
        return type == ClassificationType.AMINO_ACID_INSERTION

    def is_token_instance(self, token):
        """Return if the token is an Amino Acid Insertion token instance."""
        return isinstance(token, AminoAcidInsertionToken)
