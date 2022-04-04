"""Module for Tokenization."""
from abc import ABC, abstractmethod
from typing import Optional

from variation.schemas.token_response_schema import Token


class Tokenizer(ABC):
    """The tokenizer class."""

    @abstractmethod
    def match(self, input_string: str) -> Optional[Token]:
        """Return tokens that match the input string."""
        raise NotImplementedError
