"""A module for Protein Substitution Tokenization."""
from typing import List, Optional

from variation.schemas.token_response_schema import ProteinSubstitutionToken,\
    TokenMatchType
from .polypeptide_sequence_variation_base import PolypeptideSequenceVariationBase


class ProteinSubstitution(PolypeptideSequenceVariationBase):
    """Class for tokenizing Protein Substitution."""

    def match(self, input_string: str) -> Optional[ProteinSubstitutionToken]:
        """Return a ProteinSubstitutionToken match if one exists.

        :param str input_string: The input string to match
        :return: A ProteinSubstitutionToken if a match exists.
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

        if "." in input_string:
            if not input_string.startswith("p."):
                return None
            p_count = input_string.count("p.")
            if p_count == 1:
                psub_parts = self.splitter.split(input_string)
            elif p_count == 2:
                psub_parts = input_string.split()
        else:
            psub_parts = self.splitter.split(input_string)

        self._get_psub(psub_parts)

        if None not in self.psub.values():
            amino_acids = {self.psub["amino_acid"],
                           self.psub["new_amino_acid"]}

            if not self._is_valid_amino_acid(amino_acids):
                return None

            return ProteinSubstitutionToken(
                token=input_string,
                input_string=input_string,
                match_type=TokenMatchType.UNSPECIFIED.value,
                ref_protein=self.psub["amino_acid"],
                alt_protein=self.psub["new_amino_acid"],
                position=self.psub["position"]
            )

        return None

    def _get_psub(self, psub_parts: List) -> None:
        """Get amino acid substitution tokens.

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
