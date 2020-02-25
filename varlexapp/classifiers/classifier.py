from abc import ABC, abstractmethod
from typing import List, Optional

from ..models import Token, Classification, ClassificationType

class Classifier(ABC):

    @abstractmethod
    def match(self, tokens: List[Token]) -> Optional[Classification] :
        pass

    @abstractmethod
    def classification_type(self) -> ClassificationType:
        pass
