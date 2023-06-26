"""Module for cDNA Deletion Translation."""
from typing import Optional, List

from ga4gh.vrsatile.pydantic.vrs_models import CopyChange
from cool_seq_tool.schemas import ResidueMode

from variation.schemas.app_schemas import Endpoint
from variation.schemas.token_response_schema import AltType, CoordinateType
from variation.schemas.validation_response_schema import ValidationResult
from variation.schemas.normalize_response_schema import (
    HGVSDupDelMode as HGVSDupDelModeEnum
)
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import (
    ClassificationType, CdnaDeletionClassification
)
from variation.schemas.translation_response_schema import TranslationResult


class CdnaDeletion(Translator):
    """The cDNA Deletion Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is cDNA Deletion."""
        return type == ClassificationType.CODING_DNA_DELETION

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
        errors = []
        cds_start = validation_result.cds_start
        classification: CdnaDeletionClassification = validation_result.classification
        vrs_allele = None
        vrs_seq_loc_ac = None

        if endpoint_name == Endpoint.NORMALIZE:
            mane = await self.mane_transcript.get_mane_transcript(
                validation_result.accession, classification.pos0,
                CoordinateType.CODING_DNA,
                end_pos=classification.pos1,
                try_longest_compatible=True, residue_mode=ResidueMode.RESIDUE.value
            )

            if mane:
                vrs_seq_loc_ac = mane["refseq"]
                vrs_allele = self.vrs.to_vrs_allele(
                    vrs_seq_loc_ac, mane["pos"][0] + 1, mane["pos"][1] + 1,
                    CoordinateType.CODING_DNA, AltType.DELETION, errors,
                    cds_start=mane.get("coding_start_site", None)
                )
        else:
            vrs_seq_loc_ac = validation_result.accession
            vrs_allele = self.vrs.to_vrs_allele(
                vrs_seq_loc_ac, classification.pos0, classification.pos1,
                CoordinateType.CODING_DNA, AltType.DELETION, errors, cds_start=cds_start
            )

        if vrs_allele and vrs_seq_loc_ac:
            return TranslationResult(
                vrs_variation=vrs_allele, vrs_seq_loc_ac=vrs_seq_loc_ac
            )
        else:
            return None
