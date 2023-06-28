"""A module for Coding DNA Deletion Tokenization."""
from typing import Optional

from variation.schemas.token_response_schema import (
    CodingDNADeletionToken, CoordinateType
)
from variation.tokenizers.tokenizer import Tokenizer
from variation.regex import CNDA_GENOMIC_DELETION


class CdnaDeletion(Tokenizer):
    """Class for tokenizing Deletion at the coding dna reference sequence."""

    def match(self, input_string: str) -> Optional[CodingDNADeletionToken]:
        """Return a CodingDNADeletionToken match if one exists.

        :param input_string: The input string to match
        :return: A CodingDNADeletionToken if a match exists.
            Otherwise, None.
        """
        og_input_string = input_string
        _, input_string = self.strip_coord_prefix(
            input_string, match_coord_type=CoordinateType.CODING_DNA
        )

        if not input_string:
            return None

        match = CNDA_GENOMIC_DELETION.match(input_string)

        if match:
            match_dict = match.groupdict()

            return CodingDNADeletionToken(
                input_string=og_input_string,
                token=input_string,
                pos0=int(match_dict["pos0"]),
                pos1=int(match_dict["pos1"]) if match_dict["pos1"] else None,
                deleted_sequence=match_dict["deleted_sequence"]
            )
