"""A module for Coding DNA Silent Mutation Tokenization."""
from typing import Dict, Optional

from variation.schemas.token_response_schema import CodingDNASilentMutationToken
from .single_nucleotide_variation_base import SingleNucleotideVariationBase


class CodingDNASilentMutation(SingleNucleotideVariationBase):
    """Class for tokenizing Silent Mutation at the coding dna
    reference sequence.
    """

    def return_token(self, params: Dict) -> Optional[CodingDNASilentMutationToken]:
        """Return coding DNA silent mutation token."""
        if self.sub["coordinate_type"] == "c" and \
                self.sub["ref_nucleotide"] is None and \
                self.sub["new_nucleotide"] == "=":
            return CodingDNASilentMutationToken(**params)
