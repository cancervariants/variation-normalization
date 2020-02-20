from abc import ABC, abstractmethod

class Validator(ABC):
    @abstractmethod
    def validate(self, classification, tokens):
        pass

    @abstractmethod
    def validates_classification_type(self, classification_type):
        pass
