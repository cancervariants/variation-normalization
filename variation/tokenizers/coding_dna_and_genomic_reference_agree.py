"""A module for Reference Agree Tokenization on cDNA and genomic reference sequence."""
from typing import Optional, Union

from variation.schemas.token_response_schema import (
    CodingDNAReferenceAgreeToken, GenomicReferenceAgreeToken, CoordinateType
)
from variation.tokenizers.tokenizer import Tokenizer
from variation.regex import CDNA_GENOMIC_REFERENCE_AGREE


class CdnaGenomicReferenceAgree(Tokenizer):
    """Class for tokenizing Reference Agree on cDNA and genomic reference sequence."""

    def match(self, input_string: str) -> Optional[
        Union[CodingDNAReferenceAgreeToken, GenomicReferenceAgreeToken]
    ]:
        """Return a CodingDNAReferenceAgreeToken or GenomicReferenceAgreeToken match if
        one exists.

        :param input_string: The input string to match
        :return: A CodingDNAReferenceAgreeToken or GenomicReferenceAgreeToken if a match
            exists. Otherwise, None.
        """
        og_input_string = input_string
        coordinate_type, input_string = self.strip_coord_prefix(input_string)
        if not any((coordinate_type, input_string)):
            return None

        match = CDNA_GENOMIC_REFERENCE_AGREE.match(input_string)
        if match:
            match_dict = match.groupdict()
            params = {
                "input_string": og_input_string,
                "token": input_string,
                "coordinate_type": coordinate_type,
                "pos": int(match_dict["pos"])
            }

            if coordinate_type == CoordinateType.LINEAR_GENOMIC:
                return GenomicReferenceAgreeToken(**params)
            elif coordinate_type == CoordinateType.CODING_DNA:
                return CodingDNAReferenceAgreeToken(**params)
