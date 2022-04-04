"""A module for Genomic DelIns Tokenization."""
from typing import Dict, Optional

from variation.schemas.token_response_schema import GenomicDelInsToken
from variation.tokenizers.delins_base import DelInsBase


class GenomicDelIns(DelInsBase):
    """Class for tokenizing DelIns at the linear
    genomic reference sequence.
    """

    def return_token(self, params: Dict) -> Optional[GenomicDelInsToken]:
        """Return genomic delins token."""
        if self.parts["coordinate_type"] == "g":
            return GenomicDelInsToken(**params)
