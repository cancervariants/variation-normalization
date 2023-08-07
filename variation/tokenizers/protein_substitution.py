"""A module for Protein Substitution Tokenization."""
from typing import Optional, Union

from bioutils.sequences import aa1_to_aa3, aa3_to_aa1

from variation.regex import PROTEIN_SUBSTITUTION
from variation.schemas.token_response_schema import (
    ProteinStopGainToken,
    ProteinSubstitutionToken,
)
from variation.tokenizers.tokenizer import Tokenizer


class ProteinSubstitution(Tokenizer):
    """Class for tokenizing Protein Substitution."""

    def match(
        self, input_string: str
    ) -> Optional[Union[ProteinSubstitutionToken, ProteinStopGainToken]]:
        """Return a ProteinSubstitutionToken or ProteinStopGainToken match if one
        exists.

        :param input_string: The input string to match
        :return: A ProteinSubstitutionToken or ProteinStopGainToken if a match exists.
            Otherwise, None.
        """
        og_input_string = input_string

        if input_string.startswith(("(p.", "p.(")) and input_string.endswith(")"):
            input_string = input_string[3:-1]
        elif input_string.startswith("p."):
            input_string = input_string[2:]
        elif input_string[0] == "(" and input_string[-1] == ")":
            input_string = input_string[1:-1]

        match = PROTEIN_SUBSTITUTION.match(input_string)
        if match:
            match_dict = match.groupdict()

            ref = match_dict["ref"]
            pos = int(match_dict["pos"])
            alt = match_dict["alt"]

            # One letter codes for ref and alt
            aa1_ref = None
            aa1_alt = None

            # Ref and Alt should use the same 1 or 3 letter AA codes
            ref_upper = ref
            alt_upper = alt
            try:
                # see if it's 1 AA already
                aa1_to_aa3(ref_upper)
                aa1_to_aa3(alt_upper)
            except KeyError:
                # maybe 3 letter AA code was used
                try:
                    aa1_ref = aa3_to_aa1(ref)
                    aa1_alt = "*" if alt == "*" else aa3_to_aa1(alt)
                except KeyError:
                    pass
            else:
                aa1_ref = ref
                aa1_alt = alt

            if aa1_alt and aa1_ref:
                params = {
                    "input_string": og_input_string,
                    "token": f"{aa1_ref}{pos}{aa1_alt}",
                    "pos": pos,
                    "ref": aa1_ref,
                    "alt": aa1_alt,
                }

                if aa1_alt == "*":
                    return ProteinStopGainToken(**params)
                else:
                    return ProteinSubstitutionToken(**params)
