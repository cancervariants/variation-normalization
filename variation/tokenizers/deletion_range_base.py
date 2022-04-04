"""A module for tokenizing genomic deletion ranges."""
from typing import Dict, Optional, List
from abc import abstractmethod

from variation.schemas.token_response_schema import TokenMatchType, \
    DeletionRange
from .tokenizer import Tokenizer


class DeletionRangeBase(Tokenizer):
    """The tokenizer class for genomic deletion range."""

    def __init__(self) -> None:
        """Initialize the Genomic Deletion Range Class."""
        self.parts = None

    def match(self, input_string: str) -> Optional[DeletionRange]:
        """Return tokens that match the input string.

        :param str input_string: Input string
        :return: DeletionRange token if a match is found
        """
        if input_string is None:
            return None

        self.parts = {
            "token": input_string,
            "input_string": input_string,
            "match_type": TokenMatchType.UNSPECIFIED.value,
            "start_pos1_del": None,
            "start_pos2_del": None,
            "end_pos1_del": None,
            "end_pos2_del": None,
            "coordinate_type": None
        }

        input_string = str(input_string).lower()
        if not input_string.endswith("del"):
            return None

        if input_string.startswith("g."):
            self.parts["coordinate_type"] = "g"
        elif input_string.startswith("c."):
            self.parts["coordinate_type"] = "c"
        elif input_string.startswith("p."):
            self.parts["coordinate_type"] = "p"
        else:
            return None

        parts = input_string.split("_")
        self._get_parts(parts)
        return self.return_token(self.parts)

    @abstractmethod
    def _get_parts(self, parts: List) -> None:
        """Set `self.parts` for genomic deletion range

        :param List parts: Parts of input string
        """
        raise NotImplementedError

    @abstractmethod
    def return_token(self, params: Dict[str, str]) -> Optional[DeletionRange]:
        """Return token instance.

        :param Dict params: Params for DeletionRange token
        :return: DeletionRange token
        """
        raise NotImplementedError
