"""A module for gnomad VCF tokenization"""
from typing import Optional, List
import re
from .tokenizer import Tokenizer
from variation.schemas.token_response_schema import TokenMatchType,\
    Token, ChromosomeToken, GenomicSubstitutionToken


class GnomadVCF(Tokenizer):
    """The gnomad VCF tokenizer class"""

    def __init__(self) -> None:
        """Initialize the gnomad VCF tokenizer class"""
        self.splitter = re.compile(
            r"^(?P<chromosome>(chr)?([1-9]|[1][0-9]|[2][0-2]|X|Y))-"
            r"(?P<pos>[1-9]\d*)-(?P<ref>(?i)[actg]+)-(?P<alt>(?i)[actg]+)$")

    def match(self, input_string) -> Optional[List[Token]]:
        """Return a GnomadVCFToken if a match exists.

        :param str input_string: The input string to match
        :return: List of tokens
        """
        if input_string is None:
            return None

        match = self.splitter.match(input_string)
        if not match:
            return None

        tokens = list()
        params = match.groupdict()

        if "chr" not in params["chromosome"]:
            params["chromosome"] = f"chr{params['chromosome']}"
        else:
            params["chromosome"] = params["chromosome"].lower()
        for field in ["ref", "alt"]:
            params[field] = params[field].upper()

        tokens.append(ChromosomeToken(
            token=params["chromosome"],
            input_string=input_string,
            match_type=TokenMatchType.UNSPECIFIED.value,
            chromosome=params["chromosome"])
        )

        if len(params['ref']) == 1 and len(params['alt']) == 1:
            tokens.append(GenomicSubstitutionToken(
                token=input_string,
                input_string=input_string,
                match_type=TokenMatchType.UNSPECIFIED.value,
                position=params['pos'],
                ref_nucleotide=params['ref'],
                new_nucleotide=params['alt']
            ))

        return tokens
