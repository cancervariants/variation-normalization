from abc import ABC, abstractmethod


class Tokenizer(ABC):

    @abstractmethod
    def match(self, input_string):
        pass

