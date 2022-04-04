"""A module for Insertion Tokenization Base Class."""
from abc import abstractmethod
from typing import Optional, Dict, List

from variation.schemas.token_response_schema import Insertion, TokenMatchType, Token
from .tokenizer import Tokenizer
from .caches import AminoAcidCache, NucleotideCache
from .tokenize_base import TokenizeBase


class InsertionBase(Tokenizer):
    """Class for tokenizing Insertion."""

    def __init__(self, amino_acid_cache: AminoAcidCache,
                 nucleotide_cache: NucleotideCache) -> None:
        """Initialize the Insertion Class.

        :param AminoAcidCache amino_acid_cache: Valid amino acid codes
        :param NucleotideCache nucleotide_cache: Valid nucleotides
        """
        self.parts = None
        self.tokenize_base = TokenizeBase(amino_acid_cache, nucleotide_cache)

    def match(self, input_string: str) -> Optional[Insertion]:
        """Return tokens that match the input string."""
        if input_string is None:
            return None

        self.parts = {
            "token": input_string,
            "input_string": input_string,
            "match_type": TokenMatchType.UNSPECIFIED.value,
            "start_pos_flank": None,
            "end_pos_flank": None,
            "inserted_sequence": None,
            "inserted_sequence2": None,
            "coordinate_type": None
        }

        input_string = str(input_string).lower()

        if input_string.startswith("(") and input_string.endswith(")"):
            input_string = input_string[1:-1]

        conditions = (
            "ins" in input_string,
            "del" not in input_string and "delins" not in input_string,
            input_string.startswith("c.") or input_string.startswith("g.")
        )
        if not all(conditions):
            return None

        parts = input_string.split("ins")
        self._get_parts(parts)
        return self.return_token(self.parts)

    def _get_parts(self, parts: List) -> None:
        """Get parts for Insertion.

        :param List parts: Parts of input string
        """
        if len(parts) != 2:
            return None

        if not parts[0].startswith("c.") and not parts[0].startswith("g."):
            return None

        # Get reference sequence
        coordinate_type = parts[0][:1]
        parts[0] = parts[0][2:]

        positions = self.tokenize_base.get_positions_deleted(parts)
        if not positions:
            return None
        start_pos = positions[0]
        end_pos = positions[1]

        inserted_sequences = \
            self.tokenize_base.get_transcript_genomic_inserted_sequence(parts)
        if not inserted_sequences or not inserted_sequences[0]:
            return None

        if "_" in parts[0] and parts[0].count("_") ==\
                2 and not inserted_sequences[1]:
            return None

        self.parts["start_pos_flank"] = start_pos
        self.parts["end_pos_flank"] = end_pos
        self.parts["inserted_sequence"] = inserted_sequences[0]
        self.parts["inserted_sequence2"] = inserted_sequences[1]
        self.parts["coordinate_type"] = coordinate_type

    @abstractmethod
    def return_token(self, params: Dict[str, str]) -> Optional[Token]:
        """Return token instance."""
        raise NotImplementedError
