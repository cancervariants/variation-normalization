"""A module for DelIns Tokenization Base Class."""
import re
from abc import abstractmethod
from typing import Optional, Dict
from .tokenizer import Tokenizer
from variant.schemas.token_response_schema import DelIns, TokenMatchType
from variant.tokenizers.caches import NucleotideCache


class DelInsBase(Tokenizer):
    """Class for tokenizing DelIns."""

    def __init__(self) -> None:
        """Initialize the DelIns Base Class."""
        self.splitter = re.compile(r'delins')
        self.parts = None
        self.nucleotide_cache = NucleotideCache()

    def match(self, input_string: str) -> Optional[DelIns]:
        """Return tokens that match the input string."""
        if input_string is None:
            return None

        input_string = str(input_string).lower()

        if ('delins' not in input_string) or (
                (not input_string.startswith('c.') and (not
                 input_string.startswith('g.')))):
            return None

        self.parts = {
            'pos1_del': None,
            'pos2_del': None,
            'inserted_sequence1': None,
            'inserted_sequence2': None,
            'reference_sequence': None
        }

        parts = self.splitter.split(input_string)
        self._get_parts(parts)
        params = {
            'token': input_string,
            'input_string': input_string,
            'match_type': TokenMatchType.UNSPECIFIED.value,
            'pos1_del': self.parts['pos1_del'],
            'pos2_del': self.parts['pos2_del'],
            'inserted_sequence1': self.parts['inserted_sequence1'],
            'inserted_sequence2': self.parts['inserted_sequence2'],
            'reference_sequence': self.parts['reference_sequence']
        }
        return self.return_token(params)

    def _get_parts(self, parts):
        """Get parts for DelIns.

        :param list parts: Parts of input string
        """
        if len(parts) != 2:
            return

        if not parts[0].startswith('c.') and not parts[0].startswith('g.'):
            return

        # Get reference sequence
        reference_sequence = parts[0][:1]
        parts[0] = parts[0][2:]

        # Check positions deleted
        if '_' in parts[0] and parts[0].count('_') == 1:
            positions = parts[0].split('_')
            pos1_del = positions[0]
            pos2_del = positions[1]
            if not pos1_del.isdigit() or not pos2_del.isdigit():
                return
        else:
            pos1_del = parts[0]
            pos2_del = None
            if not pos1_del.isdigit():
                return

        # Check inserted sequences
        if '_' in parts[1] and parts[1].count('_') == 1:
            # Replaced by sequence position
            inserted_sequences = parts[1].split('_')
            inserted_sequence1 = inserted_sequences[0]
            inserted_sequence2 = inserted_sequences[1]
            if not inserted_sequence1.isdigit() or \
                    not inserted_sequence2.isdigit():
                return
        else:
            # Replaced by nucleotides
            nucleotides = ['a', 'c', 't', 'g']
            for char in parts[1]:
                if char not in nucleotides:
                    return
            inserted_sequence1 = parts[1].upper()
            inserted_sequence2 = None

        self.parts = {
            'pos1_del': pos1_del,
            'pos2_del': pos2_del,
            'inserted_sequence1': inserted_sequence1,
            'inserted_sequence2': inserted_sequence2,
            'reference_sequence': reference_sequence
        }

    @abstractmethod
    def return_token(self, params: Dict[str, str]):
        """Return token instance."""
        raise NotImplementedError
