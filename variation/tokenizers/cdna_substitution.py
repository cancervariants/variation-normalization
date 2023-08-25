"""A module for Cdna Substitution Tokenization."""
from typing import Optional

from cool_seq_tool.schemas import AnnotationLayer

from variation.regex import CDNA_GENOMIC_SUBSTITUTION
from variation.schemas.token_response_schema import CdnaSubstitutionToken
from variation.tokenizers.tokenizer import Tokenizer


class CdnaSubstitution(Tokenizer):
    """Class for tokenizing Substitution at the cdna reference sequence."""

    def match(self, input_string: str) -> Optional[CdnaSubstitutionToken]:
        """Return a CdnaSubstitutionToken match if one exists.

        :param input_string: The input string to match
        :return: A CdnaSubstitutionToken if a match exists.
            Otherwise, None.
        """
        og_input_string = input_string
        _, input_string = self.strip_coord_prefix(
            input_string, match_coord_type=AnnotationLayer.CDNA
        )
        if not input_string:
            return None

        match = CDNA_GENOMIC_SUBSTITUTION.match(input_string)

        if match:
            match_dict = match.groupdict()

            return CdnaSubstitutionToken(
                input_string=og_input_string,
                token=input_string,
                pos=int(match_dict["pos"]),
                ref=match_dict["ref"],
                alt=match_dict["alt"],
            )
