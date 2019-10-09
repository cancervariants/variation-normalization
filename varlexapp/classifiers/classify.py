from typing import Optional, List

from ..models import Token, Classification, ConfidenceRating
from . import *

class Classify:
    def __init__(self):
        self.classifiers = [
                ComplexClassifier(),
                ExpressionClassifier(),
                FusionClassifier(),
                OncogenicClassifier(),
                ProteinAlternateClassifier(),
                ProteinDelinsClassifier(),
                ProteinFrameshiftClassifier(),
                ProteinSubstitutionClassifier(),
                ProteinTerminationClassifier()
        ]


    def perform(self, tokens: List[Token]) -> List[Classification]:
        classifications = list()

        for classifier in self.classifiers:
            res = classifier.match(tokens)
            if res is not None and res.confidence == ConfidenceRating.EXACT:
                return [res]
            elif res is not None:
                classifications.append(res)

        highest_confidence = 0
        for match in classifications:
            if match.confidence.value > highest_confidence:
                highest_confidence = match.confidence.value

        return list(filter(lambda m, c=highest_confidence : m.confidence.value == c, classifications))
