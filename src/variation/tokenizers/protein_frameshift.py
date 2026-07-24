"""A module for Protein Frameshift Tokenization."""

import contextlib

from bioutils.sequences import aa1_to_aa3, aa3_to_aa1

from variation.regex import PROTEIN_FRAMESHIFT
from variation.schemas.token_response_schema import (
    ProteinFrameshiftToken,
)
from variation.tokenizers.tokenizer import Tokenizer


class ProteinFrameshift(Tokenizer):
    """Class for tokenizing Protein Frameshift."""

    def match(self, input_string: str) -> ProteinFrameshiftToken | None:
        """Return a ProteinFrameshiftToken match if one exists.

        :param input_string: The input string to match
        :return: A ProteinFrameshiftToken if a match exists. Otherwise, None.
        """
        og_input_string = input_string

        if input_string.startswith(("(p.", "p.(")) and input_string.endswith(")"):
            input_string = input_string[3:-1]
        elif input_string.startswith("p."):
            input_string = input_string[2:]
        elif input_string[0] == "(" and input_string[-1] == ")":
            input_string = input_string[1:-1]

        match = PROTEIN_FRAMESHIFT.match(input_string)
        if match:
            match_dict = match.groupdict()

            ref = match_dict["ref"]
            pos = int(match_dict["pos"])

            # One letter codes for ref
            aa1_ref = None

            # Ref should use the same 1 or 3 letter AA codes
            ref_upper = aa1_ref
            try:
                # see if it's 1 AA already
                aa1_to_aa3(ref_upper)
            except KeyError:
                # maybe 3 letter AA code was used
                with contextlib.suppress(KeyError):
                    aa1_ref = aa3_to_aa1(ref)
            else:
                aa1_ref = ref

            if aa1_ref:
                params = {
                    "input_string": og_input_string,
                    "token": f"{aa1_ref}{pos}fs",
                    "pos": pos,
                    "ref": aa1_ref,
                    "alt": "fs",
                }

                return ProteinFrameshiftToken(**params)

        return None
