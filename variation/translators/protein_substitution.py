"""Module for Protein Substitution Translation."""
from typing import Optional, List

from ga4gh.vrsatile.pydantic.vrs_models import CopyChange
from cool_seq_tool.schemas import ResidueMode

from variation.schemas.app_schemas import Endpoint
from variation.schemas.token_response_schema import AltType, CoordinateType
from variation.schemas.translation_response_schema import TranslationResult
from variation.schemas.validation_response_schema import ValidationResult
from variation.schemas.normalize_response_schema import (
    HGVSDupDelMode as HGVSDupDelModeEnum
)
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import (
    ClassificationType, ProteinSubstitutionClassification
)


class ProteinSubstitution(Translator):
    """The Protein Substitution Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Protein Substitution."""
        return type == ClassificationType.PROTEIN_SUBSTITUTION

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
        classification: ProteinSubstitutionClassification = validation_result.classification  # noqa: E501
        vrs_allele = None
        vrs_seq_loc_ac = None
        vrs_seq_loc_ac_status = "na"

        if endpoint_name == Endpoint.NORMALIZE:
            mane = await self.mane_transcript.get_mane_transcript(
                validation_result.accession, classification.pos,
                CoordinateType.PROTEIN, end_pos=classification.pos,
                try_longest_compatible=True, residue_mode=ResidueMode.RESIDUE.value
            )

            if mane:
                vrs_seq_loc_ac = mane["refseq"]
                vrs_seq_loc_ac_status = mane["status"]
                vrs_allele = self.vrs.to_vrs_allele(
                    vrs_seq_loc_ac, mane["pos"][0] + 1, mane["pos"][1] + 1,
                    CoordinateType.PROTEIN, AltType.SUBSTITUTION, warnings,
                    alt=classification.alt
                )
        else:
            vrs_seq_loc_ac = validation_result.accession
            vrs_allele = self.vrs.to_vrs_allele(
                vrs_seq_loc_ac, classification.pos, classification.pos,
                CoordinateType.PROTEIN, AltType.SUBSTITUTION, warnings,
                alt=classification.alt
            )

        if vrs_allele and vrs_seq_loc_ac:
            return TranslationResult(
                vrs_variation=vrs_allele, vrs_seq_loc_ac=vrs_seq_loc_ac,
                vrs_seq_loc_ac_status=vrs_seq_loc_ac_status,
                og_ac=validation_result.accession
            )
        else:
            return None
