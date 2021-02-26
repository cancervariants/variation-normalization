"""Module for translation."""
from variant.schemas.ga4gh_vrs import Allele
from variant.schemas.validation_response_schema import ValidationResult
from .translator import Translator
from .amino_acid_substitution import AminoAcidSubstitution
from .polypeptide_truncation import PolypeptideTruncation
from .silent_mutation import SilentMutation
from typing import List, Optional


class Translate:
    """The translation class."""

    def __init__(self) -> None:
        """Initialize the translation class."""
        self.all_translators: List[Translator] = [
            AminoAcidSubstitution(),
            PolypeptideTruncation(),
            SilentMutation()
        ]

    def perform(self, res: ValidationResult) \
            -> Optional[Allele]:
        """Translate a valid variant query."""
        for translator in self.all_translators:
            if translator.can_translate(
                    res.classification.classification_type):
                return translator.translate(res)
        return None
