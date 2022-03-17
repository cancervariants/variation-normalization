"""A module for Coding DNA Substitution Tokenization."""
from typing import Dict, Optional

from variation.schemas.token_response_schema import CodingDNASubstitutionToken
from .single_nucleotide_variation_base import SingleNucleotideVariationBase


class CodingDNASubstitution(SingleNucleotideVariationBase):
    """Class for tokenizing SNV Substitution at the coding dna
    reference sequence.
    """

    def return_token(self, params: Dict) -> Optional[CodingDNASubstitutionToken]:
        """Return coding DNA substitution token."""
        if self.sub["coordinate_type"] == "c" and \
                self.sub["ref_nucleotide"] is not None and \
                self.sub["new_nucleotide"] != "=":
            return CodingDNASubstitutionToken(**params)
