"""A module for Coding DNA Insertion Tokenization."""
from variant.schemas.token_response_schema import CodingDNAInsertionToken
from variant.tokenizers.insertion_base import InsertionBase


class CodingDNAInsertion(InsertionBase):
    """Class for tokenizing Insertion at the coding dna reference sequence."""

    def return_token(self, params):
        """Return coding DNA Insertion token."""
        if self.parts['reference_sequence'] == 'c':
            return CodingDNAInsertionToken(**params)
