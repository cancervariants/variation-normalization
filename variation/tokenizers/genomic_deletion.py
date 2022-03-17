"""A module for Genomic Deletion Tokenization."""
from typing import Dict, Optional

from variation.schemas.token_response_schema import GenomicDeletionToken
from variation.tokenizers.deletion_base import DeletionBase


class GenomicDeletion(DeletionBase):
    """Class for tokenizing Deletion at the genomic reference sequence."""

    def return_token(self, params: Dict) -> Optional[GenomicDeletionToken]:
        """Return Genomic Deletion token."""
        if self.parts["coordinate_type"] == "g":
            return GenomicDeletionToken(**params)
