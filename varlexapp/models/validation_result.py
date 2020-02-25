from typing import List

from . import Classification

class ValidationResult:
    def __init__(self,
            classification: Classification,
            is_valid: bool,
            confidence_score: float,
            concise_description: str = "",
            human_description: str = "",
            errors: List[str] = []) -> None:
        self.classification = classification
        self.is_valid = is_valid
        self.confidence_score = confidence_score
        self.concise_description = concise_description
        self.human_description = human_description
        self.errors = errors
