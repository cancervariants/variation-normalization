"""A module for Genomic Insertion Tokenization."""
from typing import Optional

from cool_seq_tool.schemas import AnnotationLayer

from variation.regex import CDNA_GENOMIC_INSERTION
from variation.schemas.token_response_schema import (
    GenomicInsertionToken,
)
from variation.tokenizers.tokenizer import Tokenizer


class GenomicInsertion(Tokenizer):
    """Class for tokenizing Insertion at the genomic reference sequence."""

    def match(self, input_string: str) -> Optional[GenomicInsertionToken]:
        """Return a GenomicInsertionToken match if one exists.

        :param input_string: The input string to match
        :return: A GenomicInsertionToken if a match exists.
            Otherwise, None.
        """
        og_input_string = input_string
        _, input_string = self.strip_coord_prefix(
            input_string, match_coord_type=AnnotationLayer.GENOMIC
        )
        if not input_string:
            return None

        match = CDNA_GENOMIC_INSERTION.match(input_string)

        if match:
            match_dict = match.groupdict()
            pos0 = int(match_dict["pos0"])
            pos1 = int(match_dict["pos1"])
            inserted_sequence = match_dict["inserted_sequence"]

            return GenomicInsertionToken(
                input_string=og_input_string,
                token=f"{pos0}_{pos1}{inserted_sequence}",
                pos0=pos0,
                pos1=pos1,
                inserted_sequence=inserted_sequence,
            )
