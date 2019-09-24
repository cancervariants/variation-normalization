from abc import ABC, abstractmethod
from typing import Optional


class Tokenizer(ABC):

    @abstractmethod
    def match(self, input_string):
        pass

