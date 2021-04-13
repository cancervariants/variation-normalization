"""A module for DNA Sequence Variant Tokenization Base Class."""
import re
from abc import abstractmethod
from typing import Optional
from .tokenizer import Tokenizer
from ..schemas.token_response_schema import Token


class DNASequenceVariantBase(Tokenizer):
    """Class for tokenizing DNA Sequence Variants."""

    def __init__(self) -> None:
        """Initialize the DNA Sequence Variant Base Class."""
        self.splitter = re.compile(r'(\d+)')
        self.sub = None

    def _set_sub(self, ref_nucleotide, position, new_nucleotide,
                 reference_sequence):
        """Initialize substitution.

        :param str ref_nucleotide: Nucleotide at reference position
        :param str position: The position nucleotide substituted
        :param str new_nucleotide: The substituted nucleotide
        :param str reference_sequence: The reference sequence used
        """
        self.sub['ref_nucleotide'] = ref_nucleotide
        self.sub['position'] = int(position)
        self.sub['new_nucleotide'] = new_nucleotide
        self.sub['reference_sequence'] = reference_sequence

    @abstractmethod
    def match(self, input_string: str) -> Optional[Token]:
        """Return tokens that match the input string."""
        raise NotImplementedError
