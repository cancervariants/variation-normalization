from abc import ABC, abstractmethod
from typing import Optional
from .token import Token


class Tokenizer(ABC):

    @abstractmethod
    def match(self, input_string):
        pass

