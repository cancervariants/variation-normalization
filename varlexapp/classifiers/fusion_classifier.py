from typing import List, Optional
from functools import reduce
import operator

from .classifier import Classifier
from ..models import Classification, Token, ConfidenceRating

class FusionClassifier(Classifier):
    def match(self, tokens: List[Token]) -> Optional[Classification]:
        token_types = list(map(lambda t: t.token_type, tokens))
        exact_matches = []
        out_of_order_matches = []
        matches_with_extra_terms = set()
        partial_matches = set()

        for candidate in self.exact_match_candidates():
            token_set = set(token_types)
            candidate_set = set(candidate)
            if token_types == candidate:
                exact_matches.append(candidate)
            elif token_set == candidate_set:
                out_of_order_matches.append(candidate)
            elif token_set > candidate_set:
                matches_with_extra_terms.update(token_set & candidate_set)
            elif len(token_set & candidate_set) > 0:
                partial_matches.update(token_set & candidate_set)

        if len(exact_matches) == 1:
            return Classification(
                    'Fusion',
                    exact_matches[0],
                    [],
                    ConfidenceRating.HIGH
            )
        elif len(out_of_order_matches) > 0:
            return Classification(
                    'Fusion',
                    list(set(reduce(operator.concat, out_of_order_matches))),
                    [],
                    ConfidenceRating.MEDIUM
            )
        elif len(matches_with_extra_terms) > 0:
            return Classification(
                    'Fusion',
                    list(matches_with_extra_terms),
                    list(token_set - matches_with_extra_terms),
                    ConfidenceRating.LOW
            )
        elif len(partial_matches) > 0:
            return Classification(
                    'Fusion',
                    list(partial_matches),
                    list(token_set - partial_matches),
                    ConfidenceRating.VERY_LOW
            )
        else:
            return None


    def exact_match_candidates(self) -> List[List[str]]:
        return [
          ['GenePair', 'Fusion'],
          ['GenePair'],
          ['GeneSymbol', 'Fusion']
        ]

