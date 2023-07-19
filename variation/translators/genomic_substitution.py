"""Module for Genomic Substitution Translation."""
from typing import Optional, List

from ga4gh.vrsatile.pydantic.vrs_models import CopyChange
from ga4gh.vrsatile.pydantic.vrsatile_models import MoleculeContext
from cool_seq_tool.schemas import ResidueMode

from variation.schemas.app_schemas import Endpoint
from variation.schemas.token_response_schema import AltType, CoordinateType
from variation.schemas.validation_response_schema import ValidationResult
from variation.schemas.normalize_response_schema import (
    HGVSDupDelMode as HGVSDupDelModeEnum
)
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import (
    ClassificationType, GenomicSubstitutionClassification,
    CdnaSubstitutionClassification
)
from variation.schemas.translation_response_schema import TranslationResult


class GenomicSubstitution(Translator):
    """The Genomic Substitution Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Genomic Substitution."""
        return type == ClassificationType.GENOMIC_SUBSTITUTION

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

        # First will translate valid result to VRS Allele
        classification: GenomicSubstitutionClassification = validation_result.classification  # noqa: E501
        vrs_allele = None
        vrs_seq_loc_ac = None
        vrs_seq_loc_ac_status = "na"

        if endpoint_name == Endpoint.NORMALIZE:
            gene = classification.gene_token.token if classification.gene_token else None  # noqa: E501
            mane = await self.mane_transcript.get_mane_transcript(
                validation_result.accession, classification.pos,
                CoordinateType.LINEAR_GENOMIC, end_pos=classification.pos,
                try_longest_compatible=True, residue_mode=ResidueMode.RESIDUE.value,
                gene=gene
            )

            if mane:
                vrs_seq_loc_ac_status = mane["status"]

                if gene:
                    if mane["strand"] == "-":
                        ref_rev = classification.ref[::-1]
                        alt_rev = classification.alt[::-1]

                        complements = {
                            "A": "T",
                            "T": "A",
                            "C": "G",
                            "G": "C"
                        }

                        ref = ""
                        alt = ""

                        for nt in ref_rev:
                            ref += complements[nt]
                        for nt in alt_rev:
                            alt += complements[nt]
                    else:
                        ref = classification.ref
                        alt = classification.alt

                    classification = CdnaSubstitutionClassification(
                        matching_tokens=classification.matching_tokens,
                        nomenclature=classification.nomenclature,
                        gene_token=classification.gene_token,
                        pos=mane["pos"][0] + 1,
                        ref=ref,
                        alt=alt,
                        so_id=classification.so_id
                    )
                    vrs_seq_loc_ac = mane["refseq"]
                    coord_type = CoordinateType.CODING_DNA
                    validation_result.classification = classification
                else:
                    vrs_seq_loc_ac = mane["alt_ac"]
                    coord_type = CoordinateType.LINEAR_GENOMIC

                vrs_allele = self.vrs.to_vrs_allele(
                    vrs_seq_loc_ac, mane["pos"][0] + 1, mane["pos"][1] + 1, coord_type,
                    AltType.SUBSTITUTION, errors, alt=classification.alt,
                    cds_start=mane["coding_start_site"] if gene else None
                )
        else:
            vrs_seq_loc_ac = validation_result.accession
            vrs_allele = self.vrs.to_vrs_allele(
                vrs_seq_loc_ac, classification.pos, classification.pos,
                CoordinateType.LINEAR_GENOMIC, AltType.SUBSTITUTION, errors,
                alt=classification.alt
            )

        if vrs_allele and vrs_seq_loc_ac:
            return TranslationResult(
                vrs_variation=vrs_allele, vrs_seq_loc_ac=vrs_seq_loc_ac,
                vrs_seq_loc_ac_status=vrs_seq_loc_ac_status,
                og_ac=validation_result.accession, validation_result=validation_result
            )
        else:
            return None
