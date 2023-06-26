"""Module for Genomic Duplication Ambiguous Translation."""
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
    AmbiguousType, ClassificationType, GenomicDuplicationAmbiguousClassification
)
from variation.schemas.translation_response_schema import TranslationResult


class GenomicDuplicationAmbiguous(Translator):
    """The Genomic Duplication Ambiguous Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Genomic Duplication Ambiguous."""
        return type == ClassificationType.GENOMIC_DUPLICATION_AMBIGUOUS

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
        classification: GenomicDuplicationAmbiguousClassification = validation_result.classification  # noqa: E501
        vrs_variation = None
        outer_coords = None

        if do_liftover or endpoint_name == Endpoint.NORMALIZE:
            errors = []
            grch38_data = await self.get_grch38_data_ambiguous(classification, errors)
            if errors:
                warnings += errors
                return None

            if classification.ambiguous_type == AmbiguousType.AMBIGUOUS_1:
                ival = self.vrs.get_ival_certain_range(
                    grch38_data["pos0"], grch38_data["pos1"], grch38_data["pos2"],
                    grch38_data["pos3"]
                )
                outer_coords = (grch38_data["pos0"], grch38_data["pos3"])
            elif classification.ambiguous_type == AmbiguousType.AMBIGUOUS_2:
                ival = models.SequenceInterval(
                    start=self.vrs.get_start_indef_range(grch38_data["pos1"]),
                    end=self.vrs.get_end_indef_range(grch38_data["pos2"]),
                    type="SequenceInterval"
                ).as_dict()
                outer_coords = (grch38_data["pos1"], grch38_data["pos2"])
            elif classification.ambiguous_type == AmbiguousType.AMBIGUOUS_5:
                ival = models.SequenceInterval(
                    start=self.vrs.get_start_indef_range(grch38_data["pos1"]),
                    end=models.Number(value=grch38_data["pos2"], type="Number"),
                    type="SequenceInterval"
                ).as_dict()
                outer_coords = (grch38_data["pos1"], grch38_data["pos2"])
            elif classification.ambiguous_type == AmbiguousType.AMBIGUOUS_7:
                ival = models.SequenceInterval(
                    start=models.Number(value=grch38_data["pos0"] - 1, type="Number"),
                    end=self.vrs.get_end_indef_range(grch38_data["pos2"]),
                    type="SequenceInterval"
                ).as_dict()
                outer_coords = (grch38_data["pos0"], grch38_data["pos2"])

            ac = grch38_data["ac"]
        else:
            ac = classification.ac

            if classification.ambiguous_type == AmbiguousType.AMBIGUOUS_1:
                ival = self.vrs.get_ival_certain_range(
                    classification.pos0, classification.pos1, classification.pos2,
                    classification.pos3
                )
                outer_coords = (classification.pos0, classification.pos3)
            elif classification.ambiguous_type == AmbiguousType.AMBIGUOUS_2:
                ival = models.SequenceInterval(
                    start=self.vrs.get_start_indef_range(classification.pos1),
                    end=self.vrs.get_end_indef_range(classification.pos2),
                    type="SequenceInterval"
                ).as_dict()
                outer_coords = (classification.pos1, classification.pos2)
            elif classification.ambiguous_type == AmbiguousType.AMBIGUOUS_5:
                ival = models.SequenceInterval(
                    start=self.vrs.get_start_indef_range(classification.pos1),
                    end=models.Number(value=classification.pos2, type="Number"),
                    type="SequenceInterval"
                ).as_dict()
                outer_coords = (classification.pos1, classification.pos2)
            elif classification.ambiguous_type == AmbiguousType.AMBIGUOUS_7:
                ival = models.SequenceInterval(
                    start=models.Number(value=classification.pos0 - 1, type="Number"),
                    end=self.vrs.get_end_indef_range(classification.pos2),
                    type="SequenceInterval"
                ).as_dict()
                outer_coords = (classification.pos0, classification.pos2)

        seq_id = self.translate_sequence_identifier(ac, errors)
        if not seq_id:
            return None

        seq_loc = self.vrs.get_sequence_loc(seq_id, ival).as_dict()

        if endpoint_name == Endpoint.NORMALIZE:
            vrs_variation = self.hgvs_dup_del_mode.interpret_variation(
                AltType.DUPLICATION_AMBIGUOUS, seq_loc, warnings, hgvs_dup_del_mode,
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
                AltType.DUPLICATION_AMBIGUOUS, outer_coords, "dup", seq_loc,
                baseline_copies=baseline_copies, copy_change=copy_change
            )

        if vrs_variation:
            return TranslationResult(vrs_variation=vrs_variation, vrs_seq_loc_ac=ac)
        else:
            return None
