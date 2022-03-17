"""A module for Genomic Silent Mutation Tokenization."""
from typing import Dict, Optional

from variation.schemas.token_response_schema import GenomicSilentMutationToken
from .single_nucleotide_variation_base import SingleNucleotideVariationBase


class GenomicSilentMutation(SingleNucleotideVariationBase):
    """Class for tokenizing Silent Mutation at the linear genomic
    reference sequence.
    """

    def return_token(self, params: Dict) -> Optional[GenomicSilentMutationToken]:
        """Return genomic silent mutation token."""
        if self.sub["coordinate_type"] == "g" and \
                self.sub["ref_nucleotide"] is None and \
                self.sub["new_nucleotide"] == "=":
            return GenomicSilentMutationToken(**params)
