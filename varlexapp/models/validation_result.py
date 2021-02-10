"""A module for the Validation Result model."""
from typing import List
from .classification import Classification
from typing import Optional


class ValidationResult:
    """The Validation Result model class."""

    def __init__(self,
                 classification: Classification,
                 is_valid: bool,
                 confidence_score: float,
                 location: Optional[dict] = None,
                 concise_description: str = "",
                 human_description: str = "",
                 errors: List[str] = []) -> None:
        """Initialize the Validation Result class."""
        self.classification = classification
        self.is_valid = is_valid
        self.confidence_score = confidence_score
        self.location = location
        self.concise_description = concise_description
        self.human_description = human_description
        self.errors = errors
