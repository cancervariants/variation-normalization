"""Module for translation."""
from typing import List, Optional

from cool_seq_tool.data_sources import MANETranscript, SeqRepoAccess, UTADatabase
from ga4gh.vrsatile.pydantic.vrs_models import CopyChange

from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.schemas.app_schemas import Endpoint
from variation.schemas.normalize_response_schema import HGVSDupDelModeOption
from variation.schemas.translation_response_schema import TranslationResult
from variation.schemas.validation_response_schema import ValidationResult
from variation.translators import (
    Amplification,
    CdnaDeletion,
    CdnaDelIns,
    CdnaInsertion,
    CdnaReferenceAgree,
    CdnaSubstitution,
    GenomicDeletion,
    GenomicDeletionAmbiguous,
    GenomicDelIns,
    GenomicDuplication,
    GenomicDuplicationAmbiguous,
    GenomicInsertion,
    GenomicReferenceAgree,
    GenomicSubstitution,
    ProteinDeletion,
    ProteinDelIns,
    ProteinInsertion,
    ProteinReferenceAgree,
    ProteinStopGain,
    ProteinSubstitution,
)
from variation.translators.translator import Translator
from variation.vrs_representation import VRSRepresentation


class Translate:
    """Class for translating to VRS representations"""

    def __init__(
        self,
        seqrepo_access: SeqRepoAccess,
        mane_transcript: MANETranscript,
        uta: UTADatabase,
        vrs: VRSRepresentation,
        hgvs_dup_del_mode: HGVSDupDelMode,
    ) -> None:
        """Initialize the Translate class. Will create an instance variable,
        `translators`, which is a list of Translator for supported variation types.


        :param seqrepo_access: Access to SeqRepo data
        :param mane_transcript: Access MANE Transcript information
        :param uta: Access to UTA queries
        :param vrs: Class for creating VRS objects
        :param hgvs_dup_del_mode: Class for interpreting HGVS duplications and deletions
        """
        params = [seqrepo_access, mane_transcript, uta, vrs, hgvs_dup_del_mode]

        self.translators: List[Translator] = [
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
            Amplification(*params),
        ]

    async def perform(
        self,
        validation_result: ValidationResult,  # this is always valid
        warnings: List[str],
        endpoint_name: Optional[Endpoint] = None,
        hgvs_dup_del_mode: HGVSDupDelModeOption = HGVSDupDelModeOption.DEFAULT,
        baseline_copies: Optional[int] = None,
        copy_change: Optional[CopyChange] = None,
        do_liftover: bool = False,
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
        translation_result = None
        for translator in self.translators:
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
                    do_liftover=do_liftover,
                )
                translation_result = result
                break

        return translation_result
