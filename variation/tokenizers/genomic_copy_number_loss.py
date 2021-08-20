"""A module for tokenizing genomic copy number loss."""
from typing import Optional
from .tokenizer import Tokenizer
from .caches import AminoAcidCache, NucleotideCache
from .tokenize_base import TokenizeBase
from variation.schemas.token_response_schema import TokenMatchType, \
    GenomicCopyNumberLossToken


class GenomicCopyNumberLoss(Tokenizer):
    """The tokenizer class for genomic copy number loss.."""

    def __init__(self, amino_acid_cache: AminoAcidCache,
                 nucleotide_cache: NucleotideCache) -> None:
        """Initialize the Genomic Copy Number Loss Class.

        :param AminoAcidCache amino_acid_cache: Valid amino acid codes
        :param NucleotideCache nucleotide_cache: Valid nucleotides
        """
        self.parts = None
        self.tokenize_base = TokenizeBase(amino_acid_cache, nucleotide_cache)

    def match(self, input_string: str) -> Optional[GenomicCopyNumberLossToken]:
        """Return tokens that match the input string."""
        if input_string is None:
            return None

        self.parts = {
            'token': input_string,
            'input_string': input_string,
            'match_type': TokenMatchType.UNSPECIFIED.value,
            'start_pos2_del': None,
            'end_pos1_del': None
        }

        input_string = str(input_string).lower()
        conditions = (
            input_string.endswith('del'),
            input_string.startswith('g.'),
            input_string.count('_') == 3
        )
        if not all(conditions):
            return None

        parts = input_string.split('_')
        self._get_parts(parts)
        if self.parts['start_pos2_del'] is None or \
                self.parts['end_pos1_del'] is None:
            return None
        return GenomicCopyNumberLossToken(**self.parts)

    def _get_parts(self, parts):
        """Set parts for genomic copy number loss.

        :param list parts: Parts of input string
        """
        conditions = (
            len(parts) == 4,
            parts[0] == 'g.(?' and parts[3] == '?)del',
            parts[1].endswith(')'),
            parts[2].startswith('(')
        )
        if all(conditions):
            parts[1] = parts[1][:-1]
            parts[2] = parts[2][1:]

            try:
                parts[1] = int(parts[1])
                parts[2] = int(parts[2])
            except ValueError:
                return None
            else:
                if parts[1] < parts[2]:
                    self.parts['start_pos2_del'] = parts[1]
                    self.parts['end_pos1_del'] = parts[2]
        return None
