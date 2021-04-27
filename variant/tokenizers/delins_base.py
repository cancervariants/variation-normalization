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

        if input_string.startswith('(') and input_string.endswith(')'):
            input_string = input_string[1:-1]

        if ('delins' not in input_string) or (
                (not input_string.startswith('c.') and (not
                 input_string.startswith('g.')))):
            return None

        self.parts = {
            'start_pos_del': None,
            'end_pos_del': None,
            'inserted_sequence1': None,
            'inserted_sequence2': None,
            'reference_sequence': None
        }

        parts = self.splitter.split(input_string)
        self._get_parts(parts)

        # TODO: implement delins range
        #  Ex: 812_829delins908_925
        if self.parts['inserted_sequence1'] is not None and \
                self.parts['inserted_sequence2'] is not None:
            return None
        params = {
            'token': input_string,
            'input_string': input_string,
            'match_type': TokenMatchType.UNSPECIFIED.value,
            'start_pos_del': self.parts['start_pos_del'],
            'end_pos_del': self.parts['end_pos_del'],
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

        positions_deleted = self._get_positions_deleted(parts)
        if not positions_deleted:
            return
        start_pos_del = positions_deleted[0]
        end_pos_del = positions_deleted[1]

        inserted_sequences = self._get_inserted_sequence(parts)
        if not inserted_sequences:
            return
        inserted_sequence1 = inserted_sequences[0]
        inserted_sequence2 = inserted_sequences[1]

        if not inserted_sequence1 and not inserted_sequence2:
            return

        self.parts = {
            'start_pos_del': start_pos_del,
            'end_pos_del': end_pos_del,
            'inserted_sequence1': inserted_sequence1,
            'inserted_sequence2': inserted_sequence2,
            'reference_sequence': reference_sequence
        }

    def _get_positions_deleted(self, parts):
        """Return position(s) deleted."""
        # Check positions deleted
        if '_' in parts[0] and parts[0].count('_') == 1:
            positions = self._get_valid_digits(parts[0])
            if not positions:
                return
            start_pos_del, end_pos_del = positions
            if start_pos_del > end_pos_del:
                return
        else:
            start_pos_del = None
            end_pos_del = parts[0]
            if not end_pos_del.isdigit():
                return None
        return start_pos_del, end_pos_del

    def _get_inserted_sequence(self, parts):
        """Return inserted sequence."""
        # Check inserted sequences
        if '_' in parts[1] and parts[1].count('_') == 1:
            # Replaced by sequence position
            inserted_sequences = self._get_valid_digits(parts[1])
            if not inserted_sequences:
                return
            inserted_sequence1, inserted_sequence2 = inserted_sequences
            if inserted_sequence1 > inserted_sequence2:
                return
        else:
            # Replaced by nucleotides
            nucleotides = ['a', 'c', 't', 'g']
            for char in parts[1]:
                if char not in nucleotides:
                    return
            inserted_sequence1 = parts[1].upper()
            inserted_sequence2 = None
        return inserted_sequence1, inserted_sequence2

    def _get_valid_digits(self, part):
        """Return valid digits after splitting on `_`."""
        digits = part.split('_')
        digit1 = digits[0]
        digit2 = digits[1]
        if not digit1.isdigit() or not digit2.isdigit():
            return None
        return digit1, digit2

    @abstractmethod
    def return_token(self, params: Dict[str, str]):
        """Return token instance."""
        raise NotImplementedError
