"""Module for Amplification Translation."""
from typing import Optional, List

from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify
from ga4gh.vrsatile.pydantic.vrs_models import CopyChange

from variation.schemas.app_schemas import Endpoint
from variation.schemas.validation_response_schema import ValidationResult
from variation.schemas.normalize_response_schema import (
    HGVSDupDelMode as HGVSDupDelModeEnum
)
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.translation_response_schema import TranslationResult
from variation.utils import get_priority_sequence_location


class Amplification(Translator):
    """The Amplification Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Amplification."""
        return type == ClassificationType.AMPLIFICATION

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
        gene_descriptor = validation_result.gene_tokens[0].gene_descriptor
        priority_seq_loc = get_priority_sequence_location(
            gene_descriptor, self.seqrepo_access
        )

        if priority_seq_loc:
            vrs_cx = models.CopyNumberChange(
                subject=models.SequenceLocation(**priority_seq_loc),
                copy_change=CopyChange.HIGH_LEVEL_GAIN.value
            )
            vrs_cx._id = ga4gh_identify(vrs_cx)
            vrs_cx = vrs_cx.as_dict()
        else:
            vrs_cx = None
            warnings.append(
                f"No VRS SequenceLocation found for gene: {gene_descriptor.gene}"
            )

        if vrs_cx:
            return TranslationResult(
                vrs_variation=vrs_cx, validation_result=validation_result
            )
        else:
            return None
