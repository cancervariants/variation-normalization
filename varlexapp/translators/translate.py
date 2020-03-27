from ..models import ValidationResult, ClassificationType, VariantRepresentation
from ..data_sources import SeqRepoAccess

from .translator import Translator
from .protein_substitution import ProteinSubstitution

from typing import List, Optional

class Translate:
    def __init__(self, seqrepo: SeqRepoAccess) -> None:
        self.seqrepo = seqrepo
        self.all_translators: List[Translator] = [
                ProteinSubstitution()
        ]

    def perform(self, res: ValidationResult) -> Optional[VariantRepresentation]:
        for translator in self.all_translators:
            if translator.can_translate(res.classification.classification_type):
                return translator.translate(res)
        return None
