"""Module for the set based classifier."""
from typing import List, Optional, Set

from variation.schemas.classification_response_schema import Classification, \
    ClassificationType, ConfidenceRating
from variation.schemas.token_response_schema import Token
from .classifier import Classifier


class SetBasedClassifier(Classifier):
    """The set based classifier class for finding classification matches."""

    def match(self, tokens: List[Token]) -> Optional[Classification]:
        """Return a classification match for a list of tokens.

        :param List[Token] tokens: List of tokens to determine classification
            match type
        :return: A classification for a list of tokens
        """
        token_types = list(map(lambda t: t.token_type, tokens))
        exact_matches: List[List[str]] = []
        out_of_order_matches: List[List[str]] = []
        matches_with_extra_terms: Set[str] = set()
        partial_matches: Set[str] = set()

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
                classification_type=self.classification_type(),
                matching_tokens=exact_matches[0],
                non_matching_tokens=[],
                all_tokens=tokens,
                confidence=ConfidenceRating.EXACT
            )
        elif len(out_of_order_matches) > 0:
            all_unordered_matches: Set[str] = set()
            for match in out_of_order_matches:
                for token in match:
                    all_unordered_matches.add(token)
            return Classification(
                classification_type=self.classification_type(),
                matching_tokens=list(all_unordered_matches),
                non_matching_tokens=[],
                all_tokens=tokens,
                confidence=ConfidenceRating.UNORDERED
            )
        elif len(matches_with_extra_terms) > 0:
            return Classification(
                classification_type=self.classification_type(),
                matching_tokens=list(matches_with_extra_terms),
                non_matching_tokens=list(token_set - matches_with_extra_terms),
                all_tokens=tokens,
                confidence=ConfidenceRating.SUPERSET
            )
        elif len(partial_matches) > 0:
            return Classification(
                classification_type=self.classification_type(),
                matching_tokens=list(partial_matches),
                non_matching_tokens=list(token_set - partial_matches),
                all_tokens=tokens,
                confidence=ConfidenceRating.INTERSECTION
            )
        else:
            return None

    def classification_type(self) -> ClassificationType:
        """Return the classification type."""
        pass

    def exact_match_candidates(self) -> List[List[str]]:
        """Return the token match candidates for a given classification.

        :return: List of tokens, where order matters, that make up a given
            classification.
        """
        pass
