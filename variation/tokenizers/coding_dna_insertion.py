"""A module for Coding DNA Insertion Tokenization."""
from typing import Dict, Optional

from variation.schemas.token_response_schema import CodingDNAInsertionToken
from variation.tokenizers.insertion_base import InsertionBase


class CodingDNAInsertion(InsertionBase):
    """Class for tokenizing Insertion at the coding dna reference sequence."""

    def return_token(self, params: Dict) -> Optional[CodingDNAInsertionToken]:
        """Return coding DNA Insertion token."""
        if self.parts["coordinate_type"] == "c":
            return CodingDNAInsertionToken(**params)
