from typing import List

from . import Classification

class ClassificationResponse:
    def __init__(self, search_term: str, classifications: List[Classification]) -> None:
        self.search_term = search_term
        self.classifications = classifications