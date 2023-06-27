"""Module for Genomic Duplication Translation."""
from typing import Optional, List

from ga4gh.vrs import models
from ga4gh.vrsatile.pydantic.vrs_models import CopyChange

from variation.schemas.app_schemas import Endpoint
from variation.schemas.token_response_schema import AltType
from variation.schemas.validation_response_schema import ValidationResult
from variation.schemas.normalize_response_schema import (
    HGVSDupDelMode as HGVSDupDelModeEnum
)
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import (
    ClassificationType, GenomicDuplicationClassification
)
from variation.schemas.translation_response_schema import TranslationResult


class GenomicDuplication(Translator):
    """The Genomic Duplication Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Genomic Duplication."""
        return type == ClassificationType.GENOMIC_DUPLICATION

    async def translate(
        self,
        validation_result: ValidationResult,
        warnings: List[str],
        endpoint_name: Optional[Endpoint] = None,
        hgvs_dup_del_mode: HGVSDupDelModeEnum = HGVSDupDelModeEnum.DEFAULT,
        baseline_copies: Optional[int] = None,
        copy_change: Optional[CopyChange] = None,
        do_liftover: bool = False
    ) -> Optional[TranslationResult]:
        """Translate to VRS Variation representation."""
        # First will translate valid result to VRS Allele
        classification: GenomicDuplicationClassification = validation_result.classification  # noqa: E501
        vrs_variation = None

        if do_liftover or endpoint_name == Endpoint.NORMALIZE:
            errors = []
            grch38_data = await self.get_grch38_data(classification, errors)
            if errors:
                warnings += errors
                return None

            pos0 = grch38_data["pos0"]
            pos1 = grch38_data["pos1"]
            ac = grch38_data["ac"]
        else:
            pos0 = classification.pos0
            pos1 = classification.pos1
            ac = validation_result.accession

        outer_coords = (pos0, pos1 if pos1 else pos0)
        ival = models.SequenceInterval(
            start=models.Number(value=pos0 - 1, type="Number"),
            end=models.Number(value=pos1 if pos1 else pos0, type="Number")
        ).as_dict()

        seq_id = self.translate_sequence_identifier(ac, warnings)
        if not seq_id:
            return None

        seq_loc = self.vrs.get_sequence_loc(seq_id, ival).as_dict()

        if endpoint_name == Endpoint.NORMALIZE:
            vrs_variation = self.hgvs_dup_del_mode.interpret_variation(
                AltType.DUPLICATION, seq_loc, warnings, hgvs_dup_del_mode, ac,
                baseline_copies=baseline_copies, copy_change=copy_change
            )
        elif endpoint_name == Endpoint.HGVS_TO_COPY_NUMBER_COUNT:
            vrs_variation = self.hgvs_dup_del_mode.copy_number_count_mode(
                "dup", seq_loc, baseline_copies
            )
        elif endpoint_name == Endpoint.HGVS_TO_COPY_NUMBER_CHANGE:
            vrs_variation = self.hgvs_dup_del_mode.copy_number_change_mode(
                "dup", seq_loc, copy_change
            )
        else:
            vrs_variation = self.hgvs_dup_del_mode.default_mode(
                AltType.DUPLICATION, outer_coords, "dup", seq_loc, ac,
                baseline_copies=baseline_copies, copy_change=copy_change
            )

        if vrs_variation:
            return TranslationResult(
                vrs_variation=vrs_variation, vrs_seq_loc_ac=ac
            )
        else:
            return None
