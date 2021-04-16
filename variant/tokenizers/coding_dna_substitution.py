"""A module for DNA Substitution Tokenization."""
from variant.schemas.token_response_schema import CodingDNASubstitutionToken
from .single_nucleotide_variant_substitution import\
    SingleNucleotideVariantSubstitution


class CodingDNASubstitution(SingleNucleotideVariantSubstitution):
    """Class for tokenizing SNV Substitution."""

    def return_token(self, params):
        """Return coding DNA substitution token."""
        if self.sub['reference_sequence'] == 'c':
            return CodingDNASubstitutionToken(**params)
