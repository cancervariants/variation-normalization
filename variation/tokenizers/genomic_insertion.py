"""A module for Genomic Insertion Tokenization."""
from typing import Dict, Optional

from variation.schemas.token_response_schema import GenomicInsertionToken
from variation.tokenizers.insertion_base import InsertionBase


class GenomicInsertion(InsertionBase):
    """Class for tokenizing Insertion at the genomic reference sequence."""

    def return_token(self, params: Dict) -> Optional[GenomicInsertionToken]:
        """Return Genomic Insertion token."""
        if self.parts["coordinate_type"] == "g":
            return GenomicInsertionToken(**params)
