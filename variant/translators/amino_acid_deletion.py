"""Module for Amino Acid Deletion Translation."""
from variant.translators.translator import Translator
from variant.schemas.classification_response_schema import ClassificationType
from variant.schemas.token_response_schema import AminoAcidDeletionToken


class AminoAcidDeletion(Translator):
    """The Amino Acid Deletion Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Amino Acid Deletion."""
        return type == ClassificationType.AMINO_ACID_DELETION

    def is_token_instance(self, token):
        """Return if the token is an Amino Acid Deletion token instance."""
        return isinstance(token, AminoAcidDeletionToken)
