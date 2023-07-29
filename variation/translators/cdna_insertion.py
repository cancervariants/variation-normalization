"""Module for Cdna insertion Translation."""
from typing import List, Optional

from cool_seq_tool.schemas import ResidueMode
from ga4gh.vrsatile.pydantic.vrs_models import CopyChange

from variation.schemas.app_schemas import Endpoint
from variation.schemas.classification_response_schema import (
    CdnaInsertionClassification,
    ClassificationType,
)
from variation.schemas.normalize_response_schema import HGVSDupDelModeOption
from variation.schemas.token_response_schema import AltType, CoordinateType
from variation.schemas.translation_response_schema import TranslationResult
from variation.schemas.validation_response_schema import ValidationResult
from variation.translators.translator import Translator


class CdnaInsertion(Translator):
    """The Cdna Insertion Translator class."""

    def can_translate(self, classification_type: ClassificationType) -> bool:
        """Determine if it's possible to translate a classification.

        :param classification_type: Classification type found
        :return: `True` if `classification_type` matches translator's classification
            type. Otherwise, `False`
        """
        return classification_type == ClassificationType.CDNA_INSERTION

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
        cds_start = validation_result.cds_start
        classification: CdnaInsertionClassification = validation_result.classification
        vrs_allele = None
        vrs_seq_loc_ac = None
        vrs_seq_loc_ac_status = "na"

        if endpoint_name == Endpoint.NORMALIZE:
            mane = await self.mane_transcript.get_mane_transcript(
                validation_result.accession,
                classification.pos0,
                CoordinateType.CDNA,
                end_pos=classification.pos1,
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
                    CoordinateType.CDNA,
                    AltType.INSERTION,
                    warnings,
                    cds_start=mane.get("coding_start_site", None),
                    alt=classification.inserted_sequence,
                )
        else:
            vrs_seq_loc_ac = validation_result.accession
            vrs_allele = self.vrs.to_vrs_allele(
                validation_result.accession,
                classification.pos0,
                classification.pos1,
                CoordinateType.CDNA,
                AltType.INSERTION,
                warnings,
                cds_start=cds_start,
                alt=classification.inserted_sequence,
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
