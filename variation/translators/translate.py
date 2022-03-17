"""Module for translation."""
from typing import List, Optional, Dict

from variation.schemas.validation_response_schema import ValidationResult
from .translator import Translator
from .protein_substitution import ProteinSubstitution
from .polypeptide_truncation import PolypeptideTruncation
from .silent_mutation import SilentMutation
from .coding_dna_substitution import CodingDNASubstitution
from .genomic_substitution import GenomicSubstitution
from .coding_dna_silent_mutation import CodingDNASilentMutation
from .genomic_silent_mutation import GenomicSilentMutation
from .protein_delins import ProteinDelIns
from .coding_dna_delins import CodingDNADelIns
from .genomic_delins import GenomicDelIns
from .protein_deletion import ProteinDeletion
from .coding_dna_deletion import CodingDNADeletion
from .genomic_deletion import GenomicDeletion
from .protein_insertion import ProteinInsertion
from .coding_dna_insertion import CodingDNAInsertion
from .genomic_insertion import GenomicInsertion
from .genomic_uncertain_deletion import GenomicUncertainDeletion
from .genomic_duplication import GenomicDuplication
from .genomic_deletion_range import GenomicDeletionRange


class Translate:
    """The translation class."""

    def __init__(self) -> None:
        """Initialize the translation class."""
        self.all_translators: List[Translator] = [
            ProteinSubstitution(),
            PolypeptideTruncation(),
            SilentMutation(),
            CodingDNASubstitution(),
            GenomicSubstitution(),
            CodingDNASilentMutation(),
            GenomicSilentMutation(),
            ProteinDelIns(),
            CodingDNADelIns(),
            GenomicDelIns(),
            ProteinDeletion(),
            CodingDNADeletion(),
            GenomicDeletion(),
            ProteinInsertion(),
            CodingDNAInsertion(),
            GenomicInsertion(),
            GenomicDeletionRange(),
            GenomicUncertainDeletion(),
            GenomicDuplication()
        ]

    def perform(self, res: ValidationResult) -> Optional[Dict]:
        """Translate a valid variation query."""
        for translator in self.all_translators:
            if translator.can_translate(
                    res.classification.classification_type):
                return translator.translate(res)
        return None
