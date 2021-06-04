"""A module for Coding DNA DelIns Tokenization."""
from variation.schemas.token_response_schema import CodingDNADelInsToken
from variation.tokenizers.delins_base import DelInsBase


class CodingDNADelIns(DelInsBase):
    """Class for tokenizing DelIns at the coding dna reference sequence."""

    def return_token(self, params):
        """Return coding DNA DelIns token."""
        if self.parts['reference_sequence'] == 'c':
            return CodingDNADelInsToken(**params)
