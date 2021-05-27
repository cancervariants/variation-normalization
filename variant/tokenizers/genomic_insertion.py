"""A module for Genomic Insertion Tokenization."""
from variant.schemas.token_response_schema import GenomicInsertionToken
from variant.tokenizers.insertion_base import InsertionBase


class GenomicInsertion(InsertionBase):
    """Class for tokenizing Insertion at the genomic reference sequence."""

    def return_token(self, params):
        """Return Genomic Insertion token."""
        if self.parts['reference_sequence'] == 'g':
            return GenomicInsertionToken(**params)
