"""Module for Genomic Duplication Translation."""
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
    ClassificationType, GenomicDuplicationClassification
)


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
    ) -> Optional[Dict]:
        """Translate to VRS Variation representation."""
        # First will translate valid result to VRS Allele
        classification: GenomicDuplicationClassification = validation_result.classification  # noqa: E501

        vrs_allele = None

        if endpoint_name == Endpoint.NORMALIZE:
            mane = await self.mane_transcript.get_mane_transcript(
                validation_result.accession, classification.pos0,
                CoordinateType.LINEAR_GENOMIC, end_pos=classification.pos1,
                try_longest_compatible=True, residue_mode=ResidueMode.RESIDUE.value
            )

            if mane:
                vrs_allele = self.vrs.to_vrs_allele(
                    mane["alt_ac"], mane["pos"][0] + 1, mane["pos"][1] + 1,
                    CoordinateType.LINEAR_GENOMIC, AltType.DUPLICATION, warnings
                )
        else:
            vrs_allele = self.vrs.to_vrs_allele(
                validation_result.accession, classification.pos0, classification.pos1,
                CoordinateType.LINEAR_GENOMIC, AltType.DUPLICATION, warnings
            )

        return vrs_allele
