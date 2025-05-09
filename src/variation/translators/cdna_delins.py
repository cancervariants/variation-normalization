"""Module for Cdna DelIns Translation."""

from cool_seq_tool.schemas import AnnotationLayer
from ga4gh.vrs import models

from variation.schemas.app_schemas import Endpoint
from variation.schemas.classification_response_schema import (
    CdnaDelInsClassification,
    ClassificationType,
)
from variation.schemas.normalize_response_schema import HGVSDupDelModeOption
from variation.schemas.token_response_schema import AltType
from variation.schemas.translation_response_schema import TranslationResult
from variation.schemas.validation_response_schema import ValidationResult
from variation.translators.translator import Translator


class CdnaDelIns(Translator):
    """The Cdna DelIns Translator class."""

    def can_translate(self, classification_type: ClassificationType) -> bool:
        """Determine if it's possible to translate a classification.

        :param classification_type: Classification type found
        :return: `True` if `classification_type` matches translator's classification
            type. Otherwise, `False`
        """
        return classification_type == ClassificationType.CDNA_DELINS

    async def translate(
        self,
        validation_result: ValidationResult,
        warnings: list[str],
        endpoint_name: Endpoint | None = None,
        hgvs_dup_del_mode: HGVSDupDelModeOption = HGVSDupDelModeOption.DEFAULT,
        baseline_copies: int | None = None,
        copy_change: models.CopyChange | None = None,
        do_liftover: bool = False,
    ) -> TranslationResult | None:
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
        cds_start = validation_result.cds_start
        classification: CdnaDelInsClassification = validation_result.classification

        return await self.get_p_or_cdna_translation_result(
            endpoint_name,
            validation_result,
            classification.pos0,
            classification.pos1,
            AltType.DELINS,
            AnnotationLayer.CDNA,
            warnings,
            cds_start=cds_start,
            alt=classification.inserted_sequence,
        )
