"""Module for Amplification Translation."""

from typing import List, Optional

from ga4gh.core import ga4gh_identify
from ga4gh.vrs import models

from variation.schemas.app_schemas import Endpoint
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.normalize_response_schema import HGVSDupDelModeOption
from variation.schemas.translation_response_schema import TranslationResult
from variation.schemas.validation_response_schema import ValidationResult
from variation.translators.translator import Translator
from variation.utils import get_priority_sequence_location


class Amplification(Translator):
    """The Amplification Translator class."""

    def can_translate(self, classification_type: ClassificationType) -> bool:
        """Determine if it's possible to translate a classification.

        :param classification_type: Classification type found
        :return: `True` if `classification_type` matches translator's classification
            type. Otherwise, `False`
        """
        return classification_type == ClassificationType.AMPLIFICATION

    async def translate(
        self,
        validation_result: ValidationResult,
        warnings: List[str],
        endpoint_name: Optional[Endpoint] = None,
        hgvs_dup_del_mode: HGVSDupDelModeOption = HGVSDupDelModeOption.DEFAULT,
        baseline_copies: Optional[int] = None,
        copy_change: Optional[models.CopyChange] = None,
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
        gene = validation_result.classification.gene_token.gene
        priority_seq_loc = get_priority_sequence_location(gene, self.seqrepo_access)

        if priority_seq_loc:
            vrs_cx = models.CopyNumberChange(
                location=models.SequenceLocation(**priority_seq_loc),
                copyChange=models.CopyChange.EFO_0030072,
            )
            vrs_cx.id = ga4gh_identify(vrs_cx)
            vrs_cx = vrs_cx.model_dump(exclude_none=True)
        else:
            vrs_cx = None
            warnings.append(f"No VRS SequenceLocation found for gene: {gene.label}")

        if vrs_cx:
            return TranslationResult(
                vrs_variation=vrs_cx, validation_result=validation_result
            )

        return None
