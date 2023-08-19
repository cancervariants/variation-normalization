"""A module for Genomic DelIns Tokenization."""
from typing import Optional

from cool_seq_tool.schemas import AnnotationLayer

from variation.regex import CDNA_GENOMIC_DELINS
from variation.schemas.token_response_schema import GenomicDelInsToken
from variation.tokenizers.tokenizer import Tokenizer


class GenomicDelIns(Tokenizer):
    """Class for tokenizing DelIns at the linear
    genomic reference sequence.
    """

    def match(self, input_string: str) -> Optional[GenomicDelInsToken]:
        """Return a GenomicDelInsToken match if one exists.

        :param input_string: The input string to match
        :return: A GenomicDelInsToken if a match exists.
            Otherwise, None.
        """
        og_input_string = input_string
        _, input_string = self.strip_coord_prefix(
            input_string, match_coord_type=AnnotationLayer.GENOMIC
        )
        if not input_string:
            return None

        match = CDNA_GENOMIC_DELINS.match(input_string)

        if match:
            match_dict = match.groupdict()

            return GenomicDelInsToken(
                input_string=og_input_string,
                token=input_string,
                pos0=int(match_dict["pos0"]),
                pos1=int(match_dict["pos1"]) if match_dict["pos1"] else None,
                inserted_sequence=match_dict["inserted_sequence"],
            )
