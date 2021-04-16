"""A module for DNA Substitution Tokenization."""
from typing import Optional, Dict
from abc import abstractmethod
from variant.schemas.token_response_schema import \
    SingleNucleotideVariantToken, TokenMatchType
from .single_nucleotide_variant_base import SingleNucleotideVariantBase


class SingleNucleotideVariantSubstitution(SingleNucleotideVariantBase):
    """Class for tokenizing SNV Substitution."""

    def match(self, input_string) -> Optional[SingleNucleotideVariantToken]:
        """Return a SingleNucleotideVariantToken match if one exists.

        :param str input_string: The input string to match
        :return: A SingleNucleotideVariantToken if a match one exists.
        """
        input_string = str(input_string).lower()
        self.sub = {
            'ref_nucleotide': None,
            'position': None,
            'new_nucleotide': None
        }

        if input_string is None:
            return None

        # TODO: Need to add m., and n.
        if 'c.' not in input_string and 'g.' not in input_string:
            return None

        sub_parts = self.splitter.split(input_string)
        self._get_sub(sub_parts)
        if None not in self.sub.values():

            params = {
                'token': input_string,
                'input_string': input_string,
                'match_type': TokenMatchType.UNSPECIFIED.value,
                'position': self.sub['position'],
                'ref_nucleotide': self.sub['ref_nucleotide'],
                'new_nucleotide': self.sub['new_nucleotide']
            }
            return self.return_token(params)
        return None

    def _get_sub(self, sub_parts):
        if len(sub_parts) == 3:
            if '(' in sub_parts[0] and ')' in sub_parts[2]:
                sub_parts[0] = sub_parts[0].split('(')[-1]
                sub_parts[2] = sub_parts[2].split(')')[0]

            if sub_parts[1].isdigit() and '>' in sub_parts[2]:
                ref_nuc, new_nuc = sub_parts[2].split('>')
                nucleotides = ['a', 't', 'c', 'g']
                if ref_nuc in nucleotides and new_nuc in nucleotides:
                    self._set_sub(ref_nuc, sub_parts[1], new_nuc,
                                  sub_parts[0].split('.')[0])

    @abstractmethod
    def return_token(self, params: Dict[str, str]):
        """Return token instance."""
        raise NotImplementedError
