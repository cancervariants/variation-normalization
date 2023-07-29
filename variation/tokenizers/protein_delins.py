"""A module for Protein DelIns Tokenization Class."""
from typing import Optional

from bioutils.sequences import aa1_to_aa3, aa3_to_aa1

from variation.regex import PROTEIN_DELINS
from variation.schemas.token_response_schema import ProteinDelInsToken
from variation.tokenizers.tokenizer import Tokenizer


class ProteinDelIns(Tokenizer):
    """Class for tokenizing DelIns on the protein reference sequence."""

    def match(self, input_string: str) -> Optional[ProteinDelInsToken]:
        """Return a ProteinDelInsToken match if one exists.

        :param input_string: The input string to match
        :return: A ProteinDelInsToken if a match exists.
            Otherwise, None.
        """
        og_input_string = input_string

        if input_string.startswith(("(p.", "p.(")) and input_string.endswith(")"):
            input_string = input_string[3:-1]
        elif input_string.startswith("p."):
            input_string = input_string[2:]
        elif input_string[0] == "(" and input_string[-1] == ")":
            input_string = input_string[1:-1]

        match = PROTEIN_DELINS.match(input_string)
        if match:
            match_dict = match.groupdict()

            aa0 = match_dict["aa0"]
            pos0 = int(match_dict["pos0"])
            aa1 = match_dict["aa1"]
            pos1 = int(match_dict["pos1"]) if match_dict["pos1"] else None
            inserted_sequence = match_dict["inserted_sequence"]

            # One letter codes for aa0, aa1, and inserted sequence
            one_letter_aa0 = None
            one_letter_aa1 = None
            one_letter_ins_seq = None

            # Should use the same 1 or 3 letter AA codes
            try:
                # see if it's 1 AA already
                aa1_to_aa3(aa0)

                if aa1:
                    aa1_to_aa3(aa1)

                if inserted_sequence:
                    aa1_to_aa3(inserted_sequence)
            except KeyError:
                # maybe 3 letter AA code was used
                try:
                    one_letter_aa0 = aa3_to_aa1(aa0)

                    if aa1:
                        one_letter_aa1 = aa3_to_aa1(aa1)

                    if inserted_sequence:
                        one_letter_ins_seq = aa3_to_aa1(inserted_sequence)
                except KeyError:
                    pass
            else:
                one_letter_aa0 = aa0
                one_letter_aa1 = aa1
                one_letter_ins_seq = inserted_sequence

            if all(
                (
                    type(aa0) == type(one_letter_aa0),
                    type(aa1) == type(one_letter_aa1),
                    type(inserted_sequence) == type(one_letter_ins_seq),
                )
            ):
                return ProteinDelInsToken(
                    input_string=og_input_string,
                    token=og_input_string,
                    aa0=one_letter_aa0,
                    pos0=pos0,
                    aa1=one_letter_aa1,
                    pos1=pos1,
                    inserted_sequence=one_letter_ins_seq,
                )
