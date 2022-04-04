"""A module for DelIns Tokenization Base Class."""
from abc import abstractmethod
from typing import Optional, Dict, List

from variation.schemas.token_response_schema import DelIns, TokenMatchType, Token
from .tokenizer import Tokenizer
from .caches import AminoAcidCache, NucleotideCache
from .tokenize_base import TokenizeBase


class DelInsBase(Tokenizer):
    """Class for tokenizing DelIns."""

    def __init__(self, amino_acid_cache: AminoAcidCache,
                 nucleotide_cache: NucleotideCache) -> None:
        """Initialize the DelIns Class.

        :param AminoAcidCache amino_acid_cache: Valid amino acid codes
        :param NucleotideCache nucleotide_cache: Valid nucleotides
        """
        self.parts = None
        self.tokenize_base = TokenizeBase(amino_acid_cache, nucleotide_cache)

    def match(self, input_string: str) -> Optional[DelIns]:
        """Return tokens that match the input string."""
        if input_string is None:
            return None

        input_string = str(input_string).lower()

        if input_string.startswith("(") and input_string.endswith(")"):
            input_string = input_string[1:-1]

        if ("delins" not in input_string) or (
                (not input_string.startswith("c.") and (not
                 input_string.startswith("g.")))):
            return None

        self.parts = {
            "start_pos_del": None,
            "end_pos_del": None,
            "inserted_sequence1": None,
            "inserted_sequence2": None,
            "coordinate_type": None
        }

        parts = input_string.split("delins")
        self._get_parts(parts)

        # TODO: implement delins range
        #  Ex: 812_829delins908_925
        if self.parts["inserted_sequence1"] is not None and \
                self.parts["inserted_sequence2"] is not None:
            return None
        params = {
            "token": input_string,
            "input_string": input_string,
            "match_type": TokenMatchType.UNSPECIFIED.value,
            "start_pos_del": self.parts["start_pos_del"],
            "end_pos_del": self.parts["end_pos_del"],
            "inserted_sequence1": self.parts["inserted_sequence1"],
            "inserted_sequence2": self.parts["inserted_sequence2"],
            "coordinate_type": self.parts["coordinate_type"]
        }
        return self.return_token(params)

    def _get_parts(self, parts: List) -> None:
        """Get parts for DelIns.

        :param List parts: Parts of input string
        """
        if len(parts) != 2:
            return

        if not parts[0].startswith("c.") and not parts[0].startswith("g."):
            return

        # Get reference sequence
        coordinate_type = parts[0][:1]
        parts[0] = parts[0][2:]

        positions_deleted = self.tokenize_base.get_positions_deleted(parts)
        if not positions_deleted:
            return
        start_pos_del = positions_deleted[0]
        end_pos_del = positions_deleted[1]

        inserted_sequences = \
            self.tokenize_base.get_transcript_genomic_inserted_sequence(parts)
        if not inserted_sequences:
            return
        inserted_sequence1 = inserted_sequences[0]
        inserted_sequence2 = inserted_sequences[1]

        if not inserted_sequence1 and not inserted_sequence2:
            return

        self.parts = {
            "start_pos_del": start_pos_del,
            "end_pos_del": end_pos_del,
            "inserted_sequence1": inserted_sequence1,
            "inserted_sequence2": inserted_sequence2,
            "coordinate_type": coordinate_type
        }

    @abstractmethod
    def return_token(self, params: Dict[str, str]) -> Optional[Token]:
        """Return token instance."""
        raise NotImplementedError
