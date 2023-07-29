"""Module for Protein Stop Gain Translation."""
from typing import Optional, List

from ga4gh.vrsatile.pydantic.vrs_models import CopyChange
from cool_seq_tool.schemas import ResidueMode

from variation.schemas.app_schemas import Endpoint
from variation.schemas.token_response_schema import AltType, CoordinateType
from variation.schemas.validation_response_schema import ValidationResult
from variation.schemas.normalize_response_schema import HGVSDupDelModeOption
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import (
    ClassificationType,
    ProteinStopGainClassification,
)
from variation.schemas.translation_response_schema import TranslationResult


class ProteinStopGain(Translator):
    """The Protein Stop Gain Translator class."""

    def can_translate(self, classification_type: ClassificationType) -> bool:
        """Determine if it's possible to translate a classification.

        :param classification_type: Classification type found
        :return: `True` if `classification_type` matches translator's classification
            type. Otherwise, `False`
        """
        return classification_type == ClassificationType.PROTEIN_STOP_GAIN

    async def translate(
        self,
        validation_result: ValidationResult,
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
        # First will translate valid result to VRS Allele
        classification: ProteinStopGainClassification = (
            validation_result.classification
        )  # noqa: E501
        vrs_allele = None
        vrs_seq_loc_ac = None
        vrs_seq_loc_ac_status = "na"

        if endpoint_name == Endpoint.NORMALIZE:
            mane = await self.mane_transcript.get_mane_transcript(
                validation_result.accession,
                classification.pos,
                CoordinateType.PROTEIN,
                end_pos=classification.pos,
                try_longest_compatible=True,
                residue_mode=ResidueMode.RESIDUE.value,
            )

            if mane:
                vrs_seq_loc_ac = mane["refseq"]
                vrs_seq_loc_ac_status = mane["status"]
                vrs_allele = self.vrs.to_vrs_allele(
                    vrs_seq_loc_ac,
                    mane["pos"][0] + 1,
                    mane["pos"][1] + 1,
                    CoordinateType.PROTEIN,
                    AltType.STOP_GAIN,
                    warnings,
                    alt=classification.alt,
                )
        else:
            vrs_seq_loc_ac = validation_result.accession
            vrs_allele = self.vrs.to_vrs_allele(
                vrs_seq_loc_ac,
                classification.pos,
                classification.pos,
                CoordinateType.PROTEIN,
                AltType.STOP_GAIN,
                warnings,
                alt=classification.alt,
            )

        if vrs_allele and vrs_seq_loc_ac:
            return TranslationResult(
                vrs_variation=vrs_allele,
                vrs_seq_loc_ac=vrs_seq_loc_ac,
                vrs_seq_loc_ac_status=vrs_seq_loc_ac_status,
                og_ac=validation_result.accession,
                validation_result=validation_result,
            )
        else:
            return None
