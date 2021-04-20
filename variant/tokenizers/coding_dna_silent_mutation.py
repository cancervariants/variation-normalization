"""A module for Coding DNA Silent Mutation Tokenization."""
from variant.schemas.token_response_schema import CodingDNASilentMutationToken
from .single_nucleotide_variant_base import SingleNucleotideVariantBase


class CodingDNASilentMutation(SingleNucleotideVariantBase):
    """Class for tokenizing SNV Substitution."""

    def return_token(self, params):
        """Return coding DNA silent mutation token."""
        if self.sub['reference_sequence'] == 'c' and \
                self.sub['ref_nucleotide'] is None:
            return CodingDNASilentMutationToken(**params)
