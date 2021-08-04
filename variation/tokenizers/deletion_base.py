"""A module for Deletion Tokenization Base Class."""
import re
from abc import abstractmethod
from typing import Optional, Dict
from .tokenizer import Tokenizer
from .caches import AminoAcidCache, NucleotideCache
from .tokenize_base import TokenizeBase
from variation.schemas.token_response_schema import Deletion, TokenMatchType


class DeletionBase(Tokenizer):
    """Class for tokenizing Deletions."""

    def __init__(self, amino_acid_cache: AminoAcidCache,
                 nucleotide_cache: NucleotideCache) -> None:
        """Initialize the Deletion Base Class.

        :param AminoAcidCache amino_acid_cache: Valid amino acid codes
        :param NucleotideCache nucleotide_cache: Valid nucleotides
        """
        self.splitter = re.compile('del')
        self.parts = None
        self.tokenize_base = TokenizeBase(amino_acid_cache, nucleotide_cache)

    def match(self, input_string: str) -> Optional[Deletion]:
        """Return tokens that match the input string."""
        if input_string is None:
            return None

        self.parts = {
            'token': input_string,
            'input_string': input_string,
            'match_type': TokenMatchType.UNSPECIFIED.value,
            'start_pos_del': None,
            'end_pos_del': None,
            'deleted_sequence': None,
            'reference_sequence': None
        }

        input_string = str(input_string).lower()
        conditions = (
            'del' in input_string,
            'ins' not in input_string and 'delins' not in input_string,
            input_string.startswith('c.') or input_string.startswith('g.')
        )
        if not all(conditions):
            return None

        parts = self.splitter.split(input_string)
        self._get_parts(parts)
        return self.return_token(self.parts)

    def _get_parts(self, parts):
        """Get parts for DelIns.

        :param list parts: Parts of input string
        """
        if len(parts) != 2:
            return

        # Get reference sequence
        reference_sequence = parts[0][:1]
        parts[0] = parts[0][2:]

        positions_deleted = self.tokenize_base.get_positions_deleted(parts)
        if not positions_deleted:
            return

        if parts[1]:
            self.parts['deleted_sequence'] = \
                self.tokenize_base.get_sequence(parts[1])

        if positions_deleted[0]:
            start_pos_del = int(positions_deleted[0])
        else:
            start_pos_del = None
        if positions_deleted[1]:
            end_pos_del = int(positions_deleted[1])
        else:
            end_pos_del = None

        self.parts['start_pos_del'] = start_pos_del
        self.parts['end_pos_del'] = end_pos_del
        self.parts['reference_sequence'] = reference_sequence

    @abstractmethod
    def return_token(self, params: Dict[str, str]):
        """Return token instance."""
        raise NotImplementedError
