"""Module for translation."""
from variant.schemas.ga4gh_vrs import Allele
from variant.schemas.validation_response_schema import ValidationResult
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
            -> Optional[Allele]:
        """Translate a valid variant query."""
        for translator in self.all_translators:
            if translator.can_translate(
                    res.classification.classification_type):
                return translator.translate(res)
        return None
