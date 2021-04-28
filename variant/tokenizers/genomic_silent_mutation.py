"""A module for Genomic Silent Mutation Tokenization."""
from variant.schemas.token_response_schema import GenomicSilentMutationToken
from .single_nucleotide_variant_base import SingleNucleotideVariantBase


class GenomicSilentMutation(SingleNucleotideVariantBase):
    """Class for tokenizing Silent Mutation at the linear genomic
    reference sequence.
    """

    def return_token(self, params):
        """Return genomic silent mutation token."""
        if self.sub['reference_sequence'] == 'g' and \
                self.sub['ref_nucleotide'] is None and \
                self.sub['new_nucleotide'] == '=':
            return GenomicSilentMutationToken(**params)
