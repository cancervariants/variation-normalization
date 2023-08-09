"""A module for Cdna Deletion Tokenization."""
from typing import Optional

from cool_seq_tool.schemas import AnnotationLayer

from variation.regex import CDNA_GENOMIC_DELINS
from variation.schemas.token_response_schema import CdnaDelInsToken
from variation.tokenizers.tokenizer import Tokenizer


class CdnaDelIns(Tokenizer):
    """Class for tokenizing delins at the cdna reference sequence."""

    def match(self, input_string: str) -> Optional[CdnaDelInsToken]:
        """Return a CdnaDelInsToken match if one exists.

        :param input_string: The input string to match
        :return: A CdnaDelInsToken if a match exists.
            Otherwise, None.
        """
        og_input_string = input_string
        _, input_string = self.strip_coord_prefix(
            input_string, match_coord_type=AnnotationLayer.CDNA
        )
        if not input_string:
            return None

        match = CDNA_GENOMIC_DELINS.match(input_string)

        if match:
            match_dict = match.groupdict()

            return CdnaDelInsToken(
                input_string=og_input_string,
                token=input_string,
                pos0=int(match_dict["pos0"]),
                pos1=int(match_dict["pos1"]) if match_dict["pos1"] else None,
                inserted_sequence=match_dict["inserted_sequence"],
            )
