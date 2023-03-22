"""A module for Deletion Tokenization Base Class."""
from abc import abstractmethod
from typing import Optional, Dict
import re

from variation.schemas.token_response_schema import Deletion, TokenMatchType
from .tokenizer import Tokenizer


class DeletionBase(Tokenizer):
    """Class for tokenizing Deletions."""

    pattern = r"^(?P<start_pos>\d+)" \
              r"(_(?P<end_pos>\d+))?del(?P<deleted_seq>[actgn]+)?$"
    splitter = re.compile(pattern)

    def match(self, input_string: str) -> Optional[Deletion]:
        """Return tokens that match the input string.

        :param str input_string: Input string
        :return: Deletion token if a match is found, else `None`
        """
        parts = {
            "token": input_string,
            "input_string": input_string,
            "match_type": TokenMatchType.UNSPECIFIED.value,
            "start_pos_del": None,
            "end_pos_del": None,
            "deleted_sequence": None,
            "coordinate_type": None
        }

        input_str_l = str(input_string).lower()

        if input_str_l.startswith("p."):
            return None

        if input_str_l.startswith(("c.", "g.")):
            parts["coordinate_type"] = input_str_l[:1]
            input_str_l = input_str_l[2:]

        if input_str_l.startswith("(") and input_str_l.endswith(")"):
            input_str_l = input_str_l[1:-1]

        match = self.splitter.match(input_str_l)
        if not match:
            return None

        params = match.groupdict()

        parts["start_pos_del"] = params["start_pos"]
        parts["end_pos_del"] = params["end_pos"]
        if parts["start_pos_del"] and parts["end_pos_del"]:
            if parts["start_pos_del"] > parts["end_pos_del"]:
                return None

        parts["deleted_sequence"] = \
            params["deleted_seq"].upper() if params["deleted_seq"] else None

        return self.return_token(parts)

    @abstractmethod
    def return_token(self, params: Dict[str, str]) -> Optional[Deletion]:
        """Return token instance.

        :param Dict params: Params for Deletion token
        :return: Deletion token
        """
        raise NotImplementedError
