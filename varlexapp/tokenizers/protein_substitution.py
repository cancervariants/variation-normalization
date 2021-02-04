import re
from typing import Optional
from .caches import AminoAcidCache
from .tokenizer import Tokenizer
from ..models import ProteinSubstitutionToken


class ProteinSubstitution(Tokenizer):
    def __init__(self, amino_acid_cache: AminoAcidCache) -> None:
        self.amino_acid_cache = amino_acid_cache
        self.splitter = re.compile(r'(\d+)')

    def match(self, input_string) -> Optional[ProteinSubstitutionToken]:
        if 'NP_' in input_string and ':p.' in input_string:
            psub_parts = self.splitter.split(input_string.split(':')[-1])
        else:
            psub_parts = self.splitter.split(input_string)

        result = self._get_psub(input_string, psub_parts)

        if result:
            return ProteinSubstitutionToken(input_string,
                                            result[0],
                                            result[1],
                                            result[2])
        return result

    def _get_psub(self, input_string, psub_parts):
        amino_acid = None
        position = None,
        new_amino_acid = None
        if len(psub_parts) == 3:
            if self._is_valid_amino_acid(psub_parts[0], psub_parts[2]):
                amino_acid = psub_parts[0]
                position = psub_parts[1]
                new_amino_acid = psub_parts[2]
            elif 'p.' in psub_parts[0]:
                # Predicted consequences given in parentheses
                if '(' in psub_parts[0] and ')' in psub_parts[-1]:
                    amino_acid = psub_parts[0].split('p.')[-1].split('(')[-1]
                    position = psub_parts[1]
                    new_amino_acid = psub_parts[-1].split(')')[0]
                else:
                    amino_acid = psub_parts[0].split('p.')[-1]
                    position = psub_parts[1]
                    new_amino_acid = psub_parts[-1]
                if not self._is_valid_amino_acid(amino_acid, new_amino_acid):
                    return None
            if amino_acid and position and new_amino_acid:
                return amino_acid, new_amino_acid, int(position)


    def _is_valid_amino_acid(self, amino_acid, new_amino_acid):
        """Returns whether or not amino acids are valid."""
        return (amino_acid in self.amino_acid_cache and
                new_amino_acid in self.amino_acid_cache)

