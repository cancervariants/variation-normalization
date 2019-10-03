from abc import ABC, abstractmethod
from typing import Set, Optional

from ..models import Token, Classification

class Classifier(ABC):

    @abstractmethod
    def match(self, tokens: Set[Token]) -> Optional[Classification] :
        pass
