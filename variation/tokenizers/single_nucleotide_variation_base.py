"""A module for Single Nucleotide Variation Tokenization Base Class."""
import re
from abc import abstractmethod
from typing import List, Optional, Dict

from variation.schemas.token_response_schema import \
    SingleNucleotideVariation, TokenMatchType
from variation.tokenizers.caches import NucleotideCache
from .tokenizer import Tokenizer


class SingleNucleotideVariationBase(Tokenizer):
    """Class for tokenizing Single Nucleotide Variations."""

    def __init__(self) -> None:
        """Initialize the Single Nucleotide Variation Base Class."""
        self.splitter = re.compile(r"(\d+)")
        self.sub = None
        self.nucleotide_cache = NucleotideCache()

    def match(self, input_string: str) -> Optional[SingleNucleotideVariation]:
        """Return a SingleNucleotideVariationToken match if one exists.

        :param str input_string: The input string to match
        :return: A SingleNucleotideVariationToken if a match one exists.
        """
        if input_string is None:
            return None

        input_string = str(input_string).lower()
        self.sub = {
            "ref_nucleotide": None,
            "position": None,
            "new_nucleotide": None,
            "coordinate_type": None
        }

        # TODO: Need to add m., and n.
        if "c." not in input_string and "g." not in input_string:
            return None

        sub_parts = self.splitter.split(input_string)
        self._get_sub(sub_parts)
        if None not in self.sub.values():
            params = {
                "token": input_string,
                "input_string": input_string,
                "match_type": TokenMatchType.UNSPECIFIED.value,
                "position": self.sub["position"],
                "ref_nucleotide": self.sub["ref_nucleotide"],
                "new_nucleotide": self.sub["new_nucleotide"]
            }
            return self.return_token(params)
        else:
            if self.sub["position"] is not None and \
                    self.sub["new_nucleotide"] == "=" and \
                    self.sub["coordinate_type"] is not None:
                params = {
                    "token": input_string,
                    "input_string": input_string,
                    "match_type": TokenMatchType.UNSPECIFIED.value,
                    "position": self.sub["position"],
                    "new_nucleotide": self.sub["new_nucleotide"]
                }
                return self.return_token(params)

        return None

    def _get_sub(self, sub_parts: List[str]) -> None:
        """Get parts for SNV.

        :param List[str] sub_parts: Split input string
        """
        if len(sub_parts) == 3:
            if "(" in sub_parts[0] and ")" in sub_parts[2]:
                sub_parts[0] = sub_parts[0].split("(")[-1]
                sub_parts[2] = sub_parts[2].split(")")[0]

            if sub_parts[1].isdigit():
                if ">" in sub_parts[2]:
                    # Substitution
                    ref_nuc, new_nuc = sub_parts[2].split(">")
                    nucleotides = self.nucleotide_cache.nucleotides.keys()
                    if ref_nuc.upper() in nucleotides \
                            and new_nuc.upper() in nucleotides:  # noqa: E501
                        self._set_sub(ref_nuc, sub_parts[1], new_nuc,
                                      sub_parts[0].split(".")[0])
                elif sub_parts[2] == "=":
                    # Silent Mutation
                    self._set_sub(None, sub_parts[1], sub_parts[2],
                                  sub_parts[0].split(".")[0])

    def _set_sub(self, ref_nucleotide: str, position: int, new_nucleotide: str,
                 coordinate_type: str) -> None:
        """Initialize substitution.

        :param str ref_nucleotide: Nucleotide at reference position
        :param str position: The position nucleotide substituted
        :param str new_nucleotide: The substituted nucleotide
        :param str coordinate_type: The coordinate type for the variation
        """
        if ref_nucleotide:
            self.sub["ref_nucleotide"] = ref_nucleotide.upper()
        self.sub["position"] = int(position)
        self.sub["new_nucleotide"] = new_nucleotide.upper()
        if coordinate_type in ["c", "g"]:
            self.sub["coordinate_type"] = coordinate_type

    @abstractmethod
    def return_token(self,
                     params: Dict[str, str]) -> Optional[SingleNucleotideVariation]:
        """Return token instance."""
        raise NotImplementedError
