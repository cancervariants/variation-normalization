"""A module for Coding DNA Deletion Tokenization."""
from variant.schemas.token_response_schema import CodingDNADeletionToken
from variant.tokenizers.deletion_base import DeletionBase


class CodingDNADeletion(DeletionBase):
    """Class for tokenizing Deletion at the coding dna reference sequence."""

    def return_token(self, params):
        """Return coding DNA Deletion token."""
        if self.parts['reference_sequence'] == 'c':
            return CodingDNADeletionToken(**params)
