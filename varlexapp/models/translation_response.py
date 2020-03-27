from . import VariantRepresentation

from typing import List

class TranslationResponse:
    def __init__(self, search_term: str, variants: List[VariantRepresentation]) -> None:
        self.search_term = search_term
        self.variants = variants

