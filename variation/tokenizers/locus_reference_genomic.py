"""Module for Locus Reference Genomic tokenization."""
import re
from typing import Optional

from variation.schemas.token_response_schema import TokenMatchType, \
    LocusReferenceGenomicToken
from .tokenizer import Tokenizer


class LocusReferenceGenomic(Tokenizer):
    """The LRG class for tokenization."""

    def __init__(self) -> None:
        """Initialize the LRG class."""
        self.splitter = re.compile(r"(\d+)")
        self.regex = r"lrg_\d+((t|p)\d+)?"
        self.parts = None

    def match(self, input_string: str) -> Optional[LocusReferenceGenomicToken]:
        """Return a LocusReferenceGenomicToken if a match exists.

        :param str input_string: The input string to match
        :return: A LocusReferenceGenomicToken if a match exists.
            Otherwise, `None`.
        """
        match = re.match(self.regex, input_string.lower())
        if not match:
            return None

        if ":" in input_string:
            input_str = input_string.lower().split(":")[0].split("lrg_")[-1]
        else:
            input_str = input_string.lower().split("lrg_")[-1]

        if input_str.isdigit():
            return LocusReferenceGenomicToken(
                token=input_string,
                input_string=input_string,
                match_type=TokenMatchType.UNSPECIFIED.value,
                id=input_str
            )

        self.parts = {
            "id": None,
            "t": None,
            "p": None,
            "token": input_string,
            "input_string": input_string,
            "match_type": TokenMatchType.UNSPECIFIED.value
        }

        valid = False
        if "t" in input_str and "p" in input_str:
            return None
        elif "t" in input_str:
            valid = self.split_t_or_p(input_str, "t")
        elif "p" in input_str:
            valid = self.split_t_or_p(input_str, "p")

        if valid:
            return LocusReferenceGenomicToken(**self.parts)

    def split_t_or_p(self, input_str: str, p_or_t: str) -> bool:
        """Split input string on `t` or `p`.

        :param str input_str: Lowered input string
        :param str p_or_t: `p` or `t`
        :return: True if valid, False otherwise
        """
        split = input_str.split(p_or_t)
        if len(split) != 2:
            return False

        self.parts["id"] = split[0]
        if not self.parts["id"].isdigit():
            return False

        self.parts[p_or_t] = split[-1]
        if not self.parts[p_or_t].isdigit():
            return False

        return True
