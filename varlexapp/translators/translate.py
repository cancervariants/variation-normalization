"""Module for translation."""
from ..models import VariantRepresentation
from varlexapp.schemas.validation_response_schema import ValidationResult
from ..data_sources import SeqRepoAccess
from .translator import Translator
from .protein_substitution import ProteinSubstitution

from typing import List, Optional


class Translate:
    """The translation class."""

    def __init__(self, seqrepo: SeqRepoAccess) -> None:
        """Initialize the translation class."""
        self.seqrepo = seqrepo
        self.all_translators: List[Translator] = [
                ProteinSubstitution()
        ]

    def perform(self, res: ValidationResult) \
            -> Optional[VariantRepresentation]:
        """Translate a valid variant query."""
        for translator in self.all_translators:
            if translator.can_translate(
                    res.classification.classification_type):
                return translator.translate(res)
        return None
