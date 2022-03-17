"""A module for Polypeptide Truncation Tokenization."""
from typing import List, Optional

from variation.schemas.token_response_schema import PolypeptideTruncationToken,\
    TokenMatchType
from .polypeptide_sequence_variation_base import PolypeptideSequenceVariationBase


class PolypeptideTruncation(PolypeptideSequenceVariationBase):
    """Class for tokenizing Protein Substitution."""

    def match(self, input_string: str) -> Optional[PolypeptideTruncationToken]:
        """Return a PolypeptideTruncationToken match if one exists.  # noqa: D202, E501

        :param str input_string: The input string to match
        :return: A PolypeptideTruncationToken if a match exists.
            Otherwise, None.
        """
        if input_string is None:
            return None

        input_string = str(input_string).lower()
        psub_parts = None
        self.psub = {
            "amino_acid": None,
            "position": None,
            "new_amino_acid": None
        }

        if input_string.startswith("(") and input_string.endswith(")"):
            input_string = input_string[1:-1]
            if input_string.endswith("*"):
                input_string = input_string.replace("*", "Ter")

        if "." in input_string:
            split_whitespace = input_string.split()
            if not input_string.startswith("p."):
                if not len(split_whitespace) == 2:
                    return None
            else:
                p_count = input_string.count("p.")
                if p_count == 1:
                    psub_parts = self.splitter.split(input_string)
            if len(split_whitespace) == 2:
                psub_parts = split_whitespace
        else:
            psub_parts = self.splitter.split(input_string)

        self._get_psub(psub_parts)

        if None not in self.psub.values():
            if self.psub["new_amino_acid"] == "*" or \
                    self.psub["new_amino_acid"] == "Ter":
                amino_acids = {self.psub["amino_acid"]}

                if not self._is_valid_amino_acid(amino_acids):
                    return None

                return PolypeptideTruncationToken(
                    token=input_string,
                    input_string=input_string,
                    match_type=TokenMatchType.UNSPECIFIED.value,
                    ref_protein=self.psub["amino_acid"],
                    position=self.psub["position"]
                )

        return None

    def _get_psub(self, psub_parts: List) -> None:
        """Get polypeptide truncation tokens.

        :param List psub_parts: The split input string
        """
        psub_parts_len = len(psub_parts)
        if psub_parts_len == 3:
            if "p." in psub_parts[0]:
                psub_parts[0] = psub_parts[0].split("p.")[-1]
            else:
                if not psub_parts[0] and psub_parts[1] and not psub_parts[2]:
                    return

            if "(" in psub_parts[0] and ")" in psub_parts[2]:
                psub_parts[0] = psub_parts[0].split("(")[-1]
                psub_parts[2] = psub_parts[2].split(")")[0]

            self._set_psub(psub_parts[0], psub_parts[1], psub_parts[2])

        elif psub_parts_len == 2:
            if "ter" in psub_parts[0]:
                if "p." not in psub_parts[0] and "p." in psub_parts[1]:
                    psub_parts[0] = f"p.{psub_parts[0]}"
                check_nonsense = f"({psub_parts[0].replace('ter', '*')})"
                if check_nonsense == psub_parts[1]:
                    psub_parts = \
                        self.splitter.split(psub_parts[0].split("p.")[-1])
                    self._get_psub(psub_parts)
