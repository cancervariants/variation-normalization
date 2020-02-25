from typing import Optional, List, Callable

from ..models import Token, Classification, ConfidenceRating
from . import *

class Classify:
    def __init__(self) -> None:
        self.classifiers: List[Classifier] = [
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
        classifications: List[Classification] = list()

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

        filter_lambda: Callable[[Classification], bool] = lambda match: match.confidence.value == highest_confidence

        return list(filter(filter_lambda, classifications))
