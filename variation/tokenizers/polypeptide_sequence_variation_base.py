"""A module for Polypeptide Sequence Variation Tokenization Base Class."""
from abc import abstractmethod
from typing import List, Optional

from bioutils.sequences import aa3_to_aa1, aa3_to_aa1_lut

from .tokenizer import Tokenizer
from ..schemas.token_response_schema import Token


class PolypeptideSequenceVariationBase(Tokenizer):
    """Class for tokenizing Polypeptide Sequence Variations."""

    valid_amino_acids = aa3_to_aa1_lut.values()

    def __init__(self) -> None:
        """Initialize the Polypeptide Sequence Variation Base Class."""
        self.psub = None

    def _set_psub(self, amino_acid: str, position: int, new_amino_acid: str) -> None:
        """Initialize protein substitution.

        :param str amino_acid: Reference amino acid
        :param int position: The position of the amino acid substituted
        :param str new_amino_acid: The new amino_acid
        """
        for (key, aa_val) in [("amino_acid", amino_acid),
                              ("new_amino_acid", new_amino_acid)]:
            if len(aa_val) == 1:
                self.psub[key] = aa_val.upper()
            else:
                try:
                    self.psub[key] = aa3_to_aa1(aa_val.capitalize())
                except KeyError:
                    self.psub[key] = None
        self.psub["position"] = int(position)

    def _is_valid_amino_acid(self, amino_acids: List) -> bool:
        """Return whether or not amino acids are valid."""
        for amino_acid_code in amino_acids:
            if amino_acid_code not in self.valid_amino_acids:
                return False
        return True

    @abstractmethod
    def match(self, input_string: str) -> Optional[Token]:
        """Return tokens that match the input string."""
        raise NotImplementedError
