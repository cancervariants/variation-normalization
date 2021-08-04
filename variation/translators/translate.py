"""Module for translation."""
from variation.schemas.ga4gh_vrs import Allele
from variation.schemas.validation_response_schema import ValidationResult
from .translator import Translator
from .amino_acid_substitution import AminoAcidSubstitution
from .polypeptide_truncation import PolypeptideTruncation
from .silent_mutation import SilentMutation
from .coding_dna_substitution import CodingDNASubstitution
from .genomic_substitution import GenomicSubstitution
from .coding_dna_silent_mutation import CodingDNASilentMutation
from .genomic_silent_mutation import GenomicSilentMutation
from .amino_acid_delins import AminoAcidDelIns
from .coding_dna_delins import CodingDNADelIns
from .genomic_delins import GenomicDelIns
from .amino_acid_deletion import AminoAcidDeletion
from .coding_dna_deletion import CodingDNADeletion
from .genomic_deletion import GenomicDeletion
from .amino_acid_insertion import AminoAcidInsertion
from .coding_dna_insertion import CodingDNAInsertion
from .genomic_insertion import GenomicInsertion
from typing import List, Optional


class Translate:
    """The translation class."""

    def __init__(self) -> None:
        """Initialize the translation class."""
        self.all_translators: List[Translator] = [
            AminoAcidSubstitution(),
            PolypeptideTruncation(),
            SilentMutation(),
            CodingDNASubstitution(),
            GenomicSubstitution(),
            CodingDNASilentMutation(),
            GenomicSilentMutation(),
            AminoAcidDelIns(),
            CodingDNADelIns(),
            GenomicDelIns(),
            AminoAcidDeletion(),
            CodingDNADeletion(),
            GenomicDeletion(),
            AminoAcidInsertion(),
            CodingDNAInsertion(),
            GenomicInsertion()
        ]

    def perform(self, res: ValidationResult) \
            -> Optional[Allele]:
        """Translate a valid variation query."""
        for translator in self.all_translators:
            if translator.can_translate(
                    res.classification.classification_type):
                return translator.translate(res)
        return None
