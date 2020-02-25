from abc import ABC, abstractmethod
from typing import Optional

from ..models import Token

class Tokenizer(ABC):

    @abstractmethod
    def match(self, input_string: str) -> Optional[Token]:
        pass

