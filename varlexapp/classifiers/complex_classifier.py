from typing import List, Optional
import operator

from .classifier import Classifier
from .set_based_classifier import SetBasedClassifier
from ..models import Classification, Token, ConfidenceRating, ClassificationType

class ComplexClassifier(SetBasedClassifier):
    def classification_type(self) -> ClassificationType:
        return ClassificationType.COMPLEX

    def exact_match_candidates(self) -> List[List[str]]:
        return [
          ['GeneSymbol', 'Amplification'],
          ['GeneSymbol', 'Exon'],
          ['GeneSymbol', 'Exon', 'Deletion'],
          ['Amplification'],
        ]

