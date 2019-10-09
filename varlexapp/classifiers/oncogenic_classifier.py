from typing import List, Optional
import operator

from .classifier import Classifier
from .set_based_classifier import SetBasedClassifier
from ..models import Classification, Token, ConfidenceRating, ClassificationType

class OncogenicClassifier(SetBasedClassifier):
    def classification_type(self) -> ClassificationType:
        return ClassificationType.ONCOGENIC

    def exact_match_candidates(self) -> List[List[str]]:
        return [
          ['GeneSymbol', 'LossOfFunction'],
          ['GeneSymbol', 'GainOfFunction'],
          ['LossOfFunction'],
          ['GainOfFunction']
        ]

