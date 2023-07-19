"""Module for Genomic Insertion Translation."""
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
    ClassificationType, GenomicInsertionClassification, CdnaInsertionClassification
)
from variation.schemas.translation_response_schema import TranslationResult


class GenomicInsertion(Translator):
    """The Genomic Insertion Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Genomic Insertion."""
        return type == ClassificationType.GENOMIC_INSERTION

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
        classification: GenomicInsertionClassification = validation_result.classification  # noqa: E501
        vrs_allele = None
        vrs_seq_loc_ac = None
        vrs_seq_loc_ac_status = "na"

        if endpoint_name == Endpoint.NORMALIZE:
            gene = classification.gene_token.token if classification.gene_token else None
            mane = await self.mane_transcript.get_mane_transcript(
                validation_result.accession, classification.pos0,
                CoordinateType.LINEAR_GENOMIC, end_pos=classification.pos1,
                try_longest_compatible=True, residue_mode=ResidueMode.RESIDUE.value,
                gene=gene
            )

            if mane:
                vrs_seq_loc_ac_status = mane["status"]
                if gene:
                    classification = CdnaInsertionClassification(
                        matching_tokens=classification.matching_tokens,
                        nomenclature=classification.nomenclature,
                        gene_token=classification.gene_token,
                        pos0=mane["pos"][0] + 1,
                        pos1=mane["pos"][1] + 1,
                        inserted_sequence=classification.inserted_sequence
                    )
                    vrs_seq_loc_ac = mane["refseq"]
                    coord_type = CoordinateType.CODING_DNA
                    validation_result.classification = classification
                else:
                    vrs_seq_loc_ac = mane["alt_ac"]
                    coord_type = CoordinateType.LINEAR_GENOMIC

                vrs_allele = self.vrs.to_vrs_allele(
                    vrs_seq_loc_ac, mane["pos"][0] + 1, mane["pos"][1] + 1, coord_type,
                    AltType.INSERTION, warnings, alt=classification.inserted_sequence,
                    cds_start=mane["coding_start_site"] if gene else None
                )
        else:
            vrs_seq_loc_ac = validation_result.accession
            vrs_allele = self.vrs.to_vrs_allele(
                vrs_seq_loc_ac, classification.pos0, classification.pos1,
                CoordinateType.LINEAR_GENOMIC, AltType.SUBSTITUTION, warnings,
                alt=classification.inserted_sequence
            )

        if vrs_allele and vrs_seq_loc_ac:
            return TranslationResult(
                vrs_variation=vrs_allele, vrs_seq_loc_ac=vrs_seq_loc_ac,
                vrs_seq_loc_ac_status=vrs_seq_loc_ac_status,
                og_ac=validation_result.accession, validation_result=validation_result
            )
        else:
            return None
