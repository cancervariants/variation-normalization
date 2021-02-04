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
        psub_parts = self.splitter.split(input_string)
        if len(psub_parts) == 3:
            amino_acid = None
            position = None,
            new_amino_acid = None
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
                return ProteinSubstitutionToken(input_string,
                                                amino_acid,
                                                new_amino_acid,
                                                int(position))
        return None

    def _is_valid_amino_acid(self, amino_acid, new_amino_acid):
        """Returns whether or not amino acids are valid."""
        return (amino_acid in self.amino_acid_cache and
                new_amino_acid in self.amino_acid_cache)

