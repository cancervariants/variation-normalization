"""A module for Genomic DelIns Tokenization."""
from variant.schemas.token_response_schema import GenomicDelInsToken
from variant.tokenizers.delins_base import DelInsBase


class GenomicDelIns(DelInsBase):
    """Class for tokenizing DelIns at the linear
    genomic reference sequence.
    """

    def return_token(self, params):
        """Return genomic delins token."""
        if self.parts['reference_sequence'] == 'g':
            return GenomicDelInsToken(**params)
