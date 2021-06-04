"""A module for Coding DNA Deletion Tokenization."""
from variation.schemas.token_response_schema import CodingDNADeletionToken
from variation.tokenizers.deletion_base import DeletionBase


class CodingDNADeletion(DeletionBase):
    """Class for tokenizing Deletion at the coding dna reference sequence."""

    def return_token(self, params):
        """Return coding DNA Deletion token."""
        if self.parts['reference_sequence'] == 'c':
            return CodingDNADeletionToken(**params)
