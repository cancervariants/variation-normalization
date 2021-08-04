"""Module for the Expression Classifier."""
from typing import List
from .set_based_classifier import SetBasedClassifier
from variation.schemas.classification_response_schema import ClassificationType


class ExpressionClassifier(SetBasedClassifier):
    """The expression classifier class."""

    def classification_type(self) -> ClassificationType:
        """Return the expression classification type."""
        return ClassificationType.EXPRESSION

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the tokens that classify as expression."""
        return [
          ['GeneSymbol', 'Expression'],
          ['GeneSymbol', 'OverExpression'],
          ['GeneSymbol', 'UnderExpression'],
          ['Expression'],
          ['OverExpression'],
          ['UnderExpression']
        ]
