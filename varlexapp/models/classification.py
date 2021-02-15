"""Module for the Classification model."""
from typing import List
from . import ClassificationType, ConfidenceRating
from varlexapp.schemas.token_response_schema import Token


class Classification:
    """The Classification class."""

    def __init__(self,
                 classification_type: ClassificationType,
                 matching_tokens: List[str],
                 non_matching_tokens: List[str],
                 all_tokens: List[Token],
                 confidence: ConfidenceRating
                 ) -> None:
        """Initialize the Classification class."""
        self.classification_type = classification_type
        self.matching_tokens = matching_tokens
        self.non_matching_tokens = non_matching_tokens
        self.all_tokens = all_tokens
        self.confidence = confidence
