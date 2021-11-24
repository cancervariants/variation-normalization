"""A module for gnomad VCF tokenization"""
from typing import Optional, List
import re
from .tokenizer import Tokenizer
from variation.schemas.token_response_schema import TokenMatchType,\
    Token, ChromosomeToken, GenomicSubstitutionToken,\
    GenomicSilentMutationToken, GenomicDeletionToken, GenomicInsertionToken


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
        params["pos"] = int(params["pos"])

        tokens.append(ChromosomeToken(
            token=params["chromosome"],
            input_string=input_string,
            match_type=TokenMatchType.UNSPECIFIED.value,
            chromosome=params["chromosome"])
        )

        ref_len = len(params["ref"])
        alt_len = len(params["alt"])

        if ref_len == alt_len:
            if ref_len == 1:
                if params["ref"] == params["alt"]:
                    tokens.append(GenomicSilentMutationToken(
                        token=input_string,
                        input_string=input_string,
                        match_type=TokenMatchType.UNSPECIFIED.value,
                        position=params["pos"],
                        ref_nucleotide=params["ref"],
                    ))
                else:
                    tokens.append(GenomicSubstitutionToken(
                        token=input_string,
                        input_string=input_string,
                        match_type=TokenMatchType.UNSPECIFIED.value,
                        position=params["pos"],
                        ref_nucleotide=params["ref"],
                        new_nucleotide=params["alt"]
                    ))
        elif ref_len < alt_len:
            # insertion
            if ref_len == 1:
                if params["ref"][0] == params["alt"][0]:
                    tokens.append(GenomicInsertionToken(
                        token=input_string,
                        input_string=input_string,
                        match_type=TokenMatchType.UNSPECIFIED.value,
                        start_pos_flank=params["pos"],
                        end_pos_flank=params["pos"] + 1,
                        inserted_sequence=params["alt"][1:]
                    ))
        else:
            # deletion
            if alt_len == 1:
                if params["ref"][0] == params["alt"][0]:
                    del_seq = params["ref"][1:]
                    tokens.append(GenomicDeletionToken(
                        token=input_string,
                        input_string=input_string,
                        match_type=TokenMatchType.UNSPECIFIED.value,
                        start_pos_del=params["pos"] + 1,
                        end_pos_del=params["pos"] + len(del_seq),
                        deleted_seqeunce=del_seq
                    ))
        return tokens
