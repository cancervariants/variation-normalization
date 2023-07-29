"""A module for Reference Agree Tokenization."""
from typing import Optional

from bioutils.sequences import aa1_to_aa3, aa3_to_aa1

from variation.regex import PROTEIN_REFERENCE_AGREE
from variation.schemas.token_response_schema import ProteinReferenceAgreeToken
from variation.tokenizers.tokenizer import Tokenizer


class ProteinReferenceAgree(Tokenizer):
    """Class for tokenizing Reference Agree on protein reference sequence."""

    def match(self, input_string: str) -> Optional[ProteinReferenceAgreeToken]:
        """Return a ProteinReferenceAgreeToken match if one exists.

        :param str input_string: The input string to match
        :return: A ProteinReferenceAgreeToken if a match exists.
            Otherwise, None.
        """
        og_input_string = input_string

        if input_string.startswith(("(p.", "p.(")) and input_string.endswith(")"):
            input_string = input_string[3:-1]
        elif input_string.startswith("p."):
            input_string = input_string[2:]
        elif input_string[0] == "(" and input_string[-1] == ")":
            input_string = input_string[1:-1]

        match = PROTEIN_REFERENCE_AGREE.match(input_string)
        if match:
            match_dict = match.groupdict()

            ref = match_dict["ref"]
            pos = int(match_dict["pos"])

            aa1_ref = None

            # Ref and Alt should use the same 1 or 3 letter AA codes
            try:
                # see if it's 1 AA already
                aa1_to_aa3(ref)
            except KeyError:
                # maybe 3 letter AA code was used
                try:
                    aa1_ref = aa3_to_aa1(ref)
                except KeyError:
                    pass
            else:
                aa1_ref = ref

            if aa1_ref:
                return ProteinReferenceAgreeToken(
                    input_string=og_input_string,
                    token=f"{aa1_ref}{pos}=",
                    pos=pos,
                    ref=aa1_ref,
                )
