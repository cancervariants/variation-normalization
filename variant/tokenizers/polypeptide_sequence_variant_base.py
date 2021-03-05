"""A module for Polypeptide Sequence Variant Tokenization Base Class."""
import re
from abc import abstractmethod
from typing import Optional
from .caches import AminoAcidCache
from .tokenizer import Tokenizer
from ..schemas.token_response_schema import Token


class PolypeptideSequenceVariantBase(Tokenizer):
    """Class for tokenizing Polypeptide Sequence Variants."""

    def __init__(self, amino_acid_cache: AminoAcidCache) -> None:
        """Initialize the Polypeptide Sequence Variant Base Class.

        :param AminoAcidCache amino_acid_cache: Valid amino acid codes.
        """
        self.amino_acid_cache = amino_acid_cache
        self.splitter = re.compile(r'(\d+)')
        self.psub = None

    def _set_psub(self, amino_acid, position, new_amino_acid):
        """Initialize protein substitution.

        :param str amino_acid: Reference amino acid
        :param str position: The position of the amino acid substituted
        :param str new_amino_acid: The new amino_acid
        """
        self.psub['amino_acid'] = amino_acid.upper() if \
            len(amino_acid) == 1 else amino_acid.capitalize()
        self.psub['position'] = int(position)
        self.psub['new_amino_acid'] = new_amino_acid if \
            len(new_amino_acid) == 1 else new_amino_acid.capitalize()

    def _is_valid_amino_acid(self, amino_acids):
        """Return whether or not amino acids are valid."""
        for amino_acid_code in amino_acids:
            if not self.amino_acid_cache.__contains__(amino_acid_code):
                return False
        return True

    @abstractmethod
    def match(self, input_string: str) -> Optional[Token]:
        """Return tokens that match the input string."""
        raise NotImplementedError
