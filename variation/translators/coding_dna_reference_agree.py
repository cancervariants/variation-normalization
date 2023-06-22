"""Module for Coding DNA Reference Agree Translation."""
from typing import Dict, Optional, List

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
    ClassificationType, CdnaReferenceAgreeClassification
)


class CdnaReferenceAgree(Translator):
    """The Coding DNA Reference Agree Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Coding DNA Reference Agree."""
        return type == ClassificationType.CODING_DNA_REFERENCE_AGREE

    async def translate(
        self,
        validation_result: ValidationResult,
        warnings: List[str],
        endpoint_name: Optional[Endpoint] = None,
        hgvs_dup_del_mode: HGVSDupDelModeEnum = HGVSDupDelModeEnum.DEFAULT,
        baseline_copies: Optional[int] = None,
        copy_change: Optional[CopyChange] = None,
        do_liftover: bool = False
    ) -> Optional[Dict]:
        """Translate to VRS Variation representation."""
        cds_start = validation_result.cds_start
        classification: CdnaReferenceAgreeClassification = validation_result.classification
        vrs_allele = None

        if endpoint_name == Endpoint.NORMALIZE:
            mane = await self.mane_transcript.get_mane_transcript(
                validation_result.accession, classification.pos,
                CoordinateType.CODING_DNA, end_pos=classification.pos,
                try_longest_compatible=True, residue_mode=ResidueMode.RESIDUE.value
            )

            vrs_allele = self.vrs.to_vrs_allele(
                mane["refseq"], mane["pos"][0] + 1, mane["pos"][1] + 1,
                CoordinateType.CODING_DNA, AltType.REFERENCE_AGREE, warnings,
                cds_start=mane.get("coding_start_site", None)
            )
        else:
            vrs_allele = self.vrs.to_vrs_allele(
                validation_result.accession, classification.pos, classification.pos,
                CoordinateType.CODING_DNA, AltType.REFERENCE_AGREE, warnings,
                cds_start=cds_start
            )

        return vrs_allele
