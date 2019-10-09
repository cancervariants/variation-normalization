from typing import List, Optional
import operator

from .classifier import Classifier
from .set_based_classifier import SetBasedClassifier
from ..models import Classification, Token, ConfidenceRating, ClassificationType

class ExpressionClassifier(SetBasedClassifier):
    def classification_type(self) -> ClassificationType:
        return ClassificationType.EXPRESSION

    def exact_match_candidates(self) -> List[List[str]]:
        return [
          ['GeneSymbol', 'Expression'],
          ['GeneSymbol', 'OverExpression'],
          ['GeneSymbol', 'UnderExpression'],
          ['Expression'],
          ['OverExpression'],
          ['UnderExpression']
        ]

