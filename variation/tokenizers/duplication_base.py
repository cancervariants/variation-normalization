"""A module for Duplication Tokenization Base Class."""
from abc import abstractmethod
from typing import Optional
from .tokenizer import Tokenizer
from variation.schemas.token_response_schema import Duplication, \
    TokenMatchType, DuplicationAltType


class DuplicationBase(Tokenizer):
    """Class for tokenizing Deletions."""

    def __init__(self) -> None:
        """Initialize the Deletion Base Class."""
        self.parts = None

    def match(self, input_string: str) -> Optional[Duplication]:
        """Return tokens that match the input string."""
        if input_string is None:
            return None

        self.parts = {
            'token': input_string,
            'input_string': input_string,
            'match_type': TokenMatchType.UNSPECIFIED.value,
            'start_pos1_dup': None,
            'start_pos2_dup': None,
            'end_pos1_dup': None,
            'end_pos2_dup': None,
            'reference_sequence': None,
            'alt_type': DuplicationAltType.DUPLICATION
        }

        input_string = str(input_string).lower()
        if not input_string.endswith('dup'):
            return None

        parts = [input_string[:-3]]
        self._get_parts(parts)
        return self.return_token()

    @abstractmethod
    def _get_parts(self, parts):
        """Get parts for DelIns.
        :param list parts: Parts of input string
        """
        raise NotImplementedError

    @abstractmethod
    def return_token(self):
        """Return token instance."""
        raise NotImplementedError