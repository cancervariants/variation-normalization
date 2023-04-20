"""A module for tokenizing Protein Deletions."""
from typing import Optional
import re

from pydantic.error_wrappers import ValidationError
from bioutils.sequences import aa3_to_aa1_lut, aa1_to_aa3

from variation.schemas.token_response_schema import ProteinDeletionToken, \
    TokenMatchType
from .tokenizer import Tokenizer


class ProteinDeletion(Tokenizer):
    """Class for tokenizing Deletions on the protein reference sequence."""

    pattern = r"^(?P<start_aa>[a-z]+)(?P<start_pos>\d+)" \
              r"(_(?P<end_aa>[a-z]+)(?P<end_pos>\d+))?del(?P<deleted_aa>[a-z]+)?$"
    splitter = re.compile(pattern)

    def match(self, input_string: str) -> Optional[ProteinDeletionToken]:
        """Return token that match the input string

        :param str input_str: Input string to tokenize
        :return: `ProteinDeletionToken` if a match is found, else `None`
        """
        parts = {
            "used_one_letter": False,
            "token": input_string,
            "input_string": input_string,
            "match_type": TokenMatchType.UNSPECIFIED.value,
            "start_aa_del": None,
            "start_pos_del": None,
            "end_aa_del": None,
            "end_pos_del": None,
            "deleted_aa": None
        }

        input_str_l = str(input_string).lower()

        if input_str_l.startswith(("c.", "g.")):
            return None

        if input_str_l.startswith("p."):
            input_str_l = input_str_l[2:]

        if input_str_l.startswith("(") and input_str_l.endswith(")"):
            input_str_l = input_str_l[1:-1]

        match = self.splitter.match(input_str_l)
        if not match:
            return None

        params = match.groupdict()
        one_letter_aa = False  # Whether or not 1 letter AA code was used

        parts["start_aa_del"] = params["start_aa"]
        parts["start_pos_del"] = params["start_pos"]
        parts["end_aa_del"] = params["end_aa"]
        parts["end_pos_del"] = params["end_pos"]
        parts["deleted_aa"] = params["deleted_aa"]

        # This ensures that start/end/deleted AA use consistent 1 or 3 AA codes
        if len(parts["start_aa_del"]) == 1:
            one_letter_aa = True

        # Validate start and end deleted amino acids (if given)
        for key in ["start_aa_del", "end_aa_del"]:
            aa = parts[key]
            if aa:
                if one_letter_aa:
                    try:
                        parts[key] = aa1_to_aa3(aa.upper())
                    except KeyError:
                        return None
                else:
                    if aa.capitalize() not in aa3_to_aa1_lut:
                        return None
                    else:
                        parts[key] = aa.capitalize()

        # Validate deleted amino acid sequence
        if parts["deleted_aa"]:
            deleted_aa = ""
            if one_letter_aa:
                for i in range(len(parts["deleted_aa"])):
                    aa = parts["deleted_aa"][i:i + 1]
                    try:
                        deleted_aa += aa3_to_aa1_lut[aa1_to_aa3(aa.upper())]
                    except KeyError:
                        return None
            else:
                for i in range(0, len(parts["deleted_aa"]), 3):
                    aa = parts["deleted_aa"][i:i + 3]
                    cap_aa = aa.capitalize()

                    if len(cap_aa) != 3 or cap_aa not in aa3_to_aa1_lut:
                        if cap_aa != "Ter":
                            return None
                    deleted_aa += cap_aa
            parts["deleted_aa"] = deleted_aa

        try:
            return ProteinDeletionToken(**parts)
        except ValidationError:
            return None
