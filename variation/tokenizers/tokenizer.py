"""Module for Tokenization."""
from abc import ABC, abstractmethod
from typing import Optional, Tuple

from cool_seq_tool.schemas import AnnotationLayer

from variation.schemas.token_response_schema import Token


class Tokenizer(ABC):
    """The tokenizer class."""

    coord_types = {k: v.value for k, v in AnnotationLayer.__members__.items()}

    @abstractmethod
    def match(self, input_string: str) -> Optional[Token]:
        """Return tokens that match the input string.

        :param input_string: Input string
        :return: Token if match was found
        """
        raise NotImplementedError

    def strip_coord_prefix(
        self, input_string: str, match_coord_type: Optional[AnnotationLayer] = None
    ) -> Tuple[Optional[AnnotationLayer], Optional[str]]:
        """Strip parentheses and coordinate type from string

        :param input_string: Input string
        :param match_coord_type: If set, the input string must have the prefix
            corresponding to this value to succeed. If this is not set, will attempt
            to find the first match of a prefix and use that as the coordinate type.
        :return: Tuple containing coordinate type for input string and stripped string,
            if successful.
        """
        coord_type = None
        stripped_str = None

        def _strip(
            coord_type: str,
            string: str,
            match_coord_type: Optional[AnnotationLayer] = None,
        ) -> str:
            """Strip parentheses and coordinate type from string

            :param input_string: Input string
            :param match_coord_type: If set, the input string must have the prefix
                corresponding to this value to succeed
            :return: Stripped string
            """
            if string.startswith(
                (f"({coord_type}.", f"{coord_type}.(")
            ) and string.endswith(")"):
                string = string[3:-1]
            elif string.startswith(f"{coord_type}."):
                string = string[2:]
            elif string[0] == "(" and string[-1] == ")":
                string = string[1:-1]
            else:
                if match_coord_type:
                    string = None

            return string

        if match_coord_type:
            coord_type = match_coord_type
            stripped_str = _strip(coord_type.value, input_string, match_coord_type)
        else:
            for k, v in self.coord_types.items():
                if f"{v}." in input_string:
                    coord_type = AnnotationLayer[k]
                    stripped_str = _strip(v, input_string)
                    break

        return coord_type, stripped_str
