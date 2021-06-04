"""A module for Genomic Silent Mutation Tokenization."""
from variation.schemas.token_response_schema import GenomicSilentMutationToken
from .single_nucleotide_variation_base import SingleNucleotideVariationBase


class GenomicSilentMutation(SingleNucleotideVariationBase):
    """Class for tokenizing Silent Mutation at the linear genomic
    reference sequence.
    """

    def return_token(self, params):
        """Return genomic silent mutation token."""
        if self.sub['reference_sequence'] == 'g' and \
                self.sub['ref_nucleotide'] is None and \
                self.sub['new_nucleotide'] == '=':
            return GenomicSilentMutationToken(**params)
