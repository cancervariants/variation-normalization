"""A module for Polypeptide Sequence Variation Tokenization Base Class."""
import re
from abc import abstractmethod
from typing import List, Optional

from .caches import AminoAcidCache
from .tokenizer import Tokenizer
from ..schemas.token_response_schema import Token


class PolypeptideSequenceVariationBase(Tokenizer):
    """Class for tokenizing Polypeptide Sequence Variations."""

    def __init__(self, amino_acid_cache: AminoAcidCache) -> None:
        """Initialize the Polypeptide Sequence Variation Base Class.

        :param AminoAcidCache amino_acid_cache: Valid amino acid codes.
        """
        self.amino_acid_cache = amino_acid_cache
        self.splitter = re.compile(r"(\d+)")
        self.psub = None

    def _set_psub(self, amino_acid: str, position: int, new_amino_acid: str) -> None:
        """Initialize protein substitution.

        :param str amino_acid: Reference amino acid
        :param int position: The position of the amino acid substituted
        :param str new_amino_acid: The new amino_acid
        """
        self.psub["amino_acid"] = amino_acid.upper() if len(amino_acid) == 1 \
            else self.amino_acid_cache.convert_three_to_one(amino_acid)
        self.psub["position"] = int(position)
        self.psub["new_amino_acid"] = new_amino_acid.upper() if \
            len(new_amino_acid) == 1 else \
            self.amino_acid_cache.convert_three_to_one(new_amino_acid)

    def _is_valid_amino_acid(self, amino_acids: List) -> bool:
        """Return whether or not amino acids are valid."""
        for amino_acid_code in amino_acids:
            if not self.amino_acid_cache.__contains__(amino_acid_code):
                return False
        return True

    @abstractmethod
    def match(self, input_string: str) -> Optional[Token]:
        """Return tokens that match the input string."""
        raise NotImplementedError
