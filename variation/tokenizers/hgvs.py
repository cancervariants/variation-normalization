"""Module for HGVS tokenization."""
import re
from typing import Optional

from cool_seq_tool.schemas import AnnotationLayer

from variation.schemas.token_response_schema import HgvsToken
from variation.tokenizers.tokenizer import Tokenizer


class HGVS(Tokenizer):
    """The HGVS tokenizer class."""

    splitter = re.compile(
        r"^(?P<accession>(NC_|NM_|NP_|ENSP|ENST)[^:\s]+):(?P<coordinate>[cgnpr])\.(?P<change>\S+)$"  # noqa: E501
    )

    def match(self, input_string: str) -> Optional[HgvsToken]:
        """Return HGVS token matches from input string.

        :param input_string: The input string to match
        :return: `HgvsToken` if HGVS match was found, else `None`
        """
        match = self.splitter.match(input_string)
        if match:
            match_dict = match.groupdict()

            return HgvsToken(
                token=input_string,
                input_string=input_string,
                accession=match_dict["accession"],
                coordinate_type=AnnotationLayer(match_dict["coordinate"]),
                change=match_dict["change"],
            )
        else:
            return None
