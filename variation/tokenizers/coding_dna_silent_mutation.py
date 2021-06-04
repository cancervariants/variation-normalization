"""A module for Coding DNA Silent Mutation Tokenization."""
from variation.schemas.token_response_schema import\
    CodingDNASilentMutationToken
from .single_nucleotide_variation_base import SingleNucleotideVariationBase


class CodingDNASilentMutation(SingleNucleotideVariationBase):
    """Class for tokenizing Silent Mutation at the coding dna
    reference sequence.
    """

    def return_token(self, params):
        """Return coding DNA silent mutation token."""
        if self.sub['reference_sequence'] == 'c' and \
                self.sub['ref_nucleotide'] is None and \
                self.sub['new_nucleotide'] == '=':
            return CodingDNASilentMutationToken(**params)
