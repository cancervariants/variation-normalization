"""A module for gnomad VCF tokenization"""
from typing import Optional, List, Dict
import re

from variation.schemas.token_response_schema import TokenMatchType, \
    Token, ChromosomeToken, GenomicSubstitutionToken, \
    GenomicSilentMutationToken, GenomicDeletionToken, GenomicInsertionToken, \
    Nomenclature
from .tokenizer import Tokenizer


class GnomadVCF(Tokenizer):
    """The gnomad VCF tokenizer class"""

    def __init__(self) -> None:
        """Initialize the gnomad VCF tokenizer class"""
        self.splitter = re.compile(
            r"^(?P<chromosome>(chr|chromosome)?([1-9]|[1][0-9]|[2][0-2]|X|Y))-"
            r"(?P<pos>[1-9]\d*)-(?P<ref>(?i)[actg]+)-(?P<alt>(?i)[actg]+)$")

    def match(self, input_string: str) -> Optional[List[Token]]:
        """Return a GnomadVCFToken if a match exists.

        :param str input_string: The input string to match
        :return: List of tokens
        """
        params = self._get_params(input_string)
        if not params:
            return None

        tokens = list()
        chromosome = params["chromosome"]
        if "chr" in chromosome:
            if "chromosome" in chromosome:
                chromosome_val = chromosome[10:]
            else:
                chromosome_val = chromosome[3:]
        else:
            chromosome_val = chromosome
        params["chromosome"] = f"chr{chromosome_val.upper()}"
        for field in ["ref", "alt"]:
            params[field] = params[field].upper()
        params["pos"] = int(params["pos"])

        tokens.append(ChromosomeToken(
            token=chromosome,
            input_string=input_string,
            match_type=TokenMatchType.UNSPECIFIED.value,
            chromosome=params["chromosome"],
            nomenclature=Nomenclature.GNOMAD_VCF)
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
                        nomenclature=Nomenclature.GNOMAD_VCF
                    ))
                else:
                    tokens.append(GenomicSubstitutionToken(
                        token=input_string,
                        input_string=input_string,
                        match_type=TokenMatchType.UNSPECIFIED.value,
                        position=params["pos"],
                        ref_nucleotide=params["ref"],
                        new_nucleotide=params["alt"],
                        nomenclature=Nomenclature.GNOMAD_VCF
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
                        inserted_sequence=params["alt"][1:],
                        nomenclature=Nomenclature.GNOMAD_VCF
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
                        deleted_sequence=del_seq,
                        nomenclature=Nomenclature.GNOMAD_VCF
                    ))
        return tokens if tokens else None

    def _get_params(self, input_string: str) -> Optional[Dict]:
        """Get chromosome, pos, ref, alt params from input string.

        :param str input_string: Input string
        :return: chromosome, pos, ref, alt data
        """
        if not input_string:
            return None
        input_string = input_string.strip()

        match = self.splitter.match(input_string)
        if not match:
            return None
        else:
            return match.groupdict()
