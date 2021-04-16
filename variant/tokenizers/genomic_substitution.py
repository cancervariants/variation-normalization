"""A module for Genomic Substitution Tokenization."""
from variant.schemas.token_response_schema import GenomicSubstitutionToken
from .single_nucleotide_variant_substitution import \
    SingleNucleotideVariantSubstitution


class GenomicSubstitution(SingleNucleotideVariantSubstitution):
    """Class for tokenizing SNV Substitution."""

    def return_token(self, params):
        """Return Genomic Substitution token."""
        if self.sub['reference_sequence'] == 'g':
            return GenomicSubstitutionToken(**params)
