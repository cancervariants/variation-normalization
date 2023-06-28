"""A module for Coding DNA Substitution Tokenization."""
from typing import Optional

from variation.schemas.token_response_schema import (
    CodingDNASubstitutionToken, CoordinateType
)
from variation.tokenizers.tokenizer import Tokenizer
from variation.regex import CDNA_GENOMIC_SUBSTITUTION


class CdnaSubstitution(Tokenizer):
    """Class for tokenizing Substitution at the coding dna reference sequence."""

    def match(self, input_string: str) -> Optional[CodingDNASubstitutionToken]:
        """Return a CodingDNASubstitutionToken match if one exists.

        :param input_string: The input string to match
        :return: A CodingDNASubstitutionToken if a match exists.
            Otherwise, None.
        """
        og_input_string = input_string
        _, input_string = self.strip_coord_prefix(
            input_string, match_coord_type=CoordinateType.CODING_DNA
        )
        if not input_string:
            return None

        match = CDNA_GENOMIC_SUBSTITUTION.match(input_string)

        if match:
            match_dict = match.groupdict()

            return CodingDNASubstitutionToken(
                input_string=og_input_string,
                token=input_string,
                pos=int(match_dict["pos"]),
                ref=match_dict["ref"],
                alt=match_dict["alt"]
            )
