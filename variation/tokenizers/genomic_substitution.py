"""A module for Genomic Substitution Tokenization."""
from typing import Dict, Optional

from variation.schemas.token_response_schema import GenomicSubstitutionToken
from .single_nucleotide_variation_base import SingleNucleotideVariationBase


class GenomicSubstitution(SingleNucleotideVariationBase):
    """Class for tokenizing SNV Substitution at the linear genomic
    reference sequence.
    """

    def return_token(self, params: Dict) -> Optional[GenomicSubstitutionToken]:
        """Return Genomic Substitution token."""
        if self.sub["coordinate_type"] == "g" and \
                self.sub["ref_nucleotide"] is not None and \
                self.sub["new_nucleotide"] != "=":
            return GenomicSubstitutionToken(**params)
