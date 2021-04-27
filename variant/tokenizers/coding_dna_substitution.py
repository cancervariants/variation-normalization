"""A module for Coding DNA Substitution Tokenization."""
from variant.schemas.token_response_schema import CodingDNASubstitutionToken
from .single_nucleotide_variant_base import SingleNucleotideVariantBase


class CodingDNASubstitution(SingleNucleotideVariantBase):
    """Class for tokenizing SNV Substitution at the coding dna
    reference sequence.
    """

    def return_token(self, params):
        """Return coding DNA substitution token."""
        if self.sub['reference_sequence'] == 'c' and \
                self.sub['ref_nucleotide'] is not None and \
                self.sub['new_nucleotide'] != '=':
            return CodingDNASubstitutionToken(**params)
