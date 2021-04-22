"""Module for Genomic Substitution Translation."""
from .dna_sequence_variant_base import DNASequenceVariantBase
from variant.schemas.classification_response_schema import ClassificationType
from variant.schemas.token_response_schema import GenomicSubstitutionToken


class GenomicSubstitution(DNASequenceVariantBase):
    """The Genomic Substitution Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Genomic Substitution."""
        return type == ClassificationType.GENOMIC_SUBSTITUTION

    def is_token_instance(self, token):
        """Return if the token is an Genomic Substitution token instance."""
        return isinstance(token, GenomicSubstitutionToken)