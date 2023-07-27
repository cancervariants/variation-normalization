"""Module for translation."""
from typing import List, Optional

from ga4gh.vrsatile.pydantic.vrs_models import CopyChange
from cool_seq_tool.data_sources import SeqRepoAccess, UTADatabase, MANETranscript

from variation.schemas.app_schemas import Endpoint
from variation.schemas.translation_response_schema import TranslationResult
from variation.schemas.validation_response_schema import ValidationResult
from variation.schemas.normalize_response_schema import HGVSDupDelModeOption
from variation.vrs_representation import VRSRepresentation
from variation.hgvs_dup_del_mode import HGVSDupDelMode
from .translator import Translator
from .protein_substitution import ProteinSubstitution
from .protein_stop_gain import ProteinStopGain
from .protein_reference_agree import ProteinReferenceAgree
from .cdna_substitution import CdnaSubstitution
from .genomic_substitution import GenomicSubstitution
from .cdna_reference_agree import CdnaReferenceAgree
from .genomic_reference_agree import GenomicReferenceAgree
from .protein_delins import ProteinDelIns
from .cdna_delins import CdnaDelIns
from .genomic_delins import GenomicDelIns
from .protein_deletion import ProteinDeletion
from .cdna_deletion import CdnaDeletion
from .genomic_deletion import GenomicDeletion
from .genomic_deletion_ambiguous import GenomicDeletionAmbiguous
from .protein_insertion import ProteinInsertion
from .cdna_insertion import CdnaInsertion
from .genomic_insertion import GenomicInsertion
from .genomic_duplication import GenomicDuplication
from .genomic_duplication_ambiguous import GenomicDuplicationAmbiguous
from .amplification import Amplification


class Translate:
    """Class for translating to VRS representations"""

    def __init__(
        self,
        seqrepo_access: SeqRepoAccess,
        mane_transcript: MANETranscript,
        uta: UTADatabase,
        vrs: VRSRepresentation,
        hgvs_dup_del_mode: HGVSDupDelMode
    ) -> None:
        """Initialize the Translate class.

        :param seqrepo_access: Access to SeqRepo data
        :param mane_transcript: Access MANE Transcript information
        :param uta: Access to UTA queries
        :param vrs: Class for creating VRS objects
        :param hgvs_dup_del_mode: Class for interpreting HGVS duplications and deletions
        """
        params = [
            seqrepo_access, mane_transcript, uta, vrs, hgvs_dup_del_mode
        ]

        self.all_translators: List[Translator] = [
            ProteinSubstitution(*params),
            CdnaSubstitution(*params),
            GenomicSubstitution(*params),
            ProteinStopGain(*params),
            ProteinReferenceAgree(*params),
            CdnaReferenceAgree(*params),
            GenomicReferenceAgree(*params),
            ProteinDelIns(*params),
            CdnaDelIns(*params),
            GenomicDelIns(*params),
            ProteinDeletion(*params),
            CdnaDeletion(*params),
            GenomicDeletion(*params),
            GenomicDeletionAmbiguous(*params),
            ProteinInsertion(*params),
            CdnaInsertion(*params),
            GenomicInsertion(*params),
            GenomicDuplication(*params),
            GenomicDuplicationAmbiguous(*params),
            Amplification(*params)
        ]

    async def perform(
        self,
        validation_result: ValidationResult,  # this is always valid
        warnings: List[str],
        endpoint_name: Optional[Endpoint] = None,
        hgvs_dup_del_mode: HGVSDupDelModeOption = HGVSDupDelModeOption.DEFAULT,
        baseline_copies: Optional[int] = None,
        copy_change: Optional[CopyChange] = None,
        do_liftover: bool = False
    ) -> Optional[TranslationResult]:
        """Translate validation result to VRS representation

        :param validation_result: Validation result for a classification
        :param endpoint_name: Name of endpoint that is being used
        :param hgvs_dup_del_mode: Mode to use for interpreting HGVS duplications and
            deletions
        :param baseline_copies: The baseline copies for a copy number count variation
        :param copy_change: The change for a copy number change variation
        :param do_liftover: Whether or not to liftover to GRCh38 assembly
        :return: Translation result if translation was successful. If translation was
            not successful, `None`
        """
        for translator in self.all_translators:
            if translator.can_translate(
                validation_result.classification.classification_type
            ):
                result = await translator.translate(
                    validation_result,
                    warnings,
                    endpoint_name=endpoint_name,
                    hgvs_dup_del_mode=hgvs_dup_del_mode,
                    baseline_copies=baseline_copies,
                    copy_change=copy_change,
                    do_liftover=do_liftover
                )
                return result
        return None
