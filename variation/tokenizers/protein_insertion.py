"""A module for Protein Insertion Tokenization Class."""
from typing import List, Optional

from pydantic.error_wrappers import ValidationError

from variation.schemas.token_response_schema import ProteinInsertionToken, \
    TokenMatchType
from .caches import AminoAcidCache, NucleotideCache
from .tokenizer import Tokenizer
from .tokenize_base import TokenizeBase


class ProteinInsertion(Tokenizer):
    """Class for tokenizing Insertions on the protein reference sequence."""

    def __init__(self, amino_acid_cache: AminoAcidCache,
                 nucleotide_cache: NucleotideCache) -> None:
        """Initialize the Protein Insertion Class.

        :param AminoAcidCache amino_acid_cache: Valid amino acid codes
        :param NucleotideCache nucleotide_cache: Valid nucleotides
        """
        self.parts = None
        self.tokenize_base = TokenizeBase(amino_acid_cache, nucleotide_cache)

    def match(self, input_string: str) -> Optional[ProteinInsertionToken]:
        """Return token that match the input string."""
        if input_string is None:
            return None

        self.parts = {
            "used_one_letter": False,
            "token": input_string,
            "input_string": input_string,
            "match_type": TokenMatchType.UNSPECIFIED.value,
            "start_aa_flank": None,
            "start_pos_flank": None,
            "end_aa_flank": None,
            "end_pos_flank": None
        }

        input_string = str(input_string).lower()

        if "c." in input_string or "g." in input_string:
            return None

        if input_string.startswith("p."):
            input_string = input_string[2:]

        if input_string.startswith("(") and input_string.endswith(")"):
            input_string = input_string[1:-1]

        if "ins" not in input_string:
            return None

        parts = input_string.split("ins")
        self._get_parts(parts)

        try:
            return ProteinInsertionToken(**self.parts)
        except ValidationError:
            return None

    def _get_parts(self, parts: List) -> None:
        """Get parts for Protein Insertion.

        :param List parts: Parts of input string
        """
        if len(parts) != 2:
            return

        range_aa_pos = self.tokenize_base.get_aa_pos_range(parts)
        if range_aa_pos:
            self.parts["start_aa_flank"] = range_aa_pos[0]
            self.parts["end_aa_flank"] = range_aa_pos[1]
            self.parts["start_pos_flank"] = range_aa_pos[2]
            self.parts["end_pos_flank"] = range_aa_pos[3]
            self.parts["used_one_letter"] = range_aa_pos[4]
        self.parts["inserted_sequence"] = \
            self.tokenize_base.get_protein_inserted_sequence(
                parts, self.parts["used_one_letter"])
