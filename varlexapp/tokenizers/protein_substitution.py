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
        if (len(psub_parts) == 3 and
                psub_parts[0] in self.amino_acid_cache and
                psub_parts[2] in self.amino_acid_cache):
            return ProteinSubstitutionToken(input_string, psub_parts[0], psub_parts[2], int(psub_parts[1]))
        else:
            return None

