"""A module for gnomad VCF tokenization"""
import re
from typing import Optional

from variation.schemas.token_response_schema import GnomadVcfToken
from variation.tokenizers.tokenizer import Tokenizer


class GnomadVCF(Tokenizer):
    """The gnomad VCF tokenizer class"""

    splitter = re.compile(
        r"^(chr|chromosome)?(?P<chromosome>([1-9]|[1][0-9]|[2][0-2]|X|Y))-"
        r"(?P<pos>[1-9]\d*)-(?P<ref>[actg]+)-(?P<alt>[actg]+)$",
        re.IGNORECASE,
    )

    def match(self, input_string: str) -> Optional[GnomadVcfToken]:
        """Return a GnomadVCFToken if a match exists.

        :param input_string: The input string to match
        :return: `Token` if gnomAD VCF match was found, else `None`
        """
        match = self.splitter.match(input_string)
        if match:
            match_dict = match.groupdict()
            chromosome = match_dict["chromosome"].upper()
            pos = int(match_dict["pos"])
            ref = match_dict["ref"].upper()
            alt = match_dict["alt"].upper()

            return GnomadVcfToken(
                token=f"{chromosome}-{pos}-{ref}-{alt}",
                input_string=input_string,
                chromosome=chromosome,
                pos=pos,
                ref=ref,
                alt=alt,
            )

        return None
