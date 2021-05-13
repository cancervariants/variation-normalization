"""A module for Genomic Deletion Tokenization."""
from variant.schemas.token_response_schema import GenomicDeletionToken
from variant.tokenizers.deletion_base import DeletionBase


class GenomicDeletion(DeletionBase):
    """Class for tokenizing Deletion at the genomic reference sequence."""

    def return_token(self, params):
        """Return Genomic Deletion token."""
        if self.parts['reference_sequence'] == 'g':
            return GenomicDeletionToken(**params)
