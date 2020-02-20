from typing import List, Optional
from functools import reduce
import operator

from .classifier import Classifier
from ..models import Classification, Token, ConfidenceRating, ClassificationType

class SetBasedClassifier(Classifier):
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
                    self.classification_type(),
                    exact_matches[0],
                    [],
                    tokens,
                    ConfidenceRating.EXACT
            )
        elif len(out_of_order_matches) > 0:
            return Classification(
                    self.classification_type(),
                    list(set(reduce(operator.concat, out_of_order_matches))),
                    [],
                    tokens,
                    ConfidenceRating.UNORDERED
            )
        elif len(matches_with_extra_terms) > 0:
            return Classification(
                    self.classification_type(),
                    list(matches_with_extra_terms),
                    list(token_set - matches_with_extra_terms),
                    tokens,
                    ConfidenceRating.SUPERSET
            )
        elif len(partial_matches) > 0:
            return Classification(
                    self.classification_type(),
                    list(partial_matches),
                    list(token_set - partial_matches),
                    tokens,
                    ConfidenceRating.INTERSECTION
            )
        else:
            return None

    def classification_type(self) -> ClassificationType:
        pass

    def exact_match_candidates(self) -> List[List[str]]:
        pass

