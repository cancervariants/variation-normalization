"""Module for Genomic Duplication Translation."""
from typing import Optional, List

from ga4gh.vrs import models
from ga4gh.vrsatile.pydantic.vrs_models import CopyChange
from ga4gh.vrsatile.pydantic.vrsatile_models import MoleculeContext
from cool_seq_tool.schemas import ResidueMode

from variation.schemas.app_schemas import Endpoint
from variation.schemas.service_schema import ClinVarAssembly
from variation.schemas.token_response_schema import AltType
from variation.schemas.validation_response_schema import ValidationResult
from variation.schemas.normalize_response_schema import (
    HGVSDupDelMode as HGVSDupDelModeEnum
)
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import (
    ClassificationType, GenomicDuplicationClassification
)
from variation.schemas.translation_response_schema import TranslationResult
from variation.utils import get_assembly


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
    ) -> Optional[TranslationResult]:
        """Translate to VRS Variation representation."""
        # First will translate valid result to VRS Allele
        classification: GenomicDuplicationClassification = validation_result.classification  # noqa: E501
        vrs_variation = None

        # TODO: Clean this up
        if do_liftover or endpoint_name == Endpoint.NORMALIZE:
            errors = []

            # Check if we need to liftover
            assembly, w = get_assembly(self.seqrepo_access, validation_result.accession)
            if w:
                warnings.append(w)
                return None
            else:
                # assembly is either 37 or 38
                if assembly == ClinVarAssembly.GRCH37:
                    grch38_data = await self.get_grch38_data(
                        classification, errors, validation_result.accession
                    )
                    if errors:
                        warnings += errors
                        return None

                    pos0 = grch38_data["pos0"]
                    pos1 = grch38_data["pos1"]
                    ac = grch38_data["ac"]
                else:
                    pos0 = classification.pos0
                    pos1 = classification.pos1
                    ac = validation_result.accession

                assembly = ClinVarAssembly.GRCH38
        else:
            pos0 = classification.pos0
            pos1 = classification.pos1
            ac = validation_result.accession
            assembly = None

        if classification.gene:
            errors = []
            if not assembly:
                grch38_data = await self.get_grch38_data(
                    classification, errors, ac
                )
                if errors:
                    warnings += errors
                    return None

                self.is_valid(
                    classification.gene, grch38_data["ac"], grch38_data["pos0"],
                    grch38_data["pos1"], errors
                )
            else:
                self.is_valid(classification.gene, ac, pos0, pos1, errors)

            if errors:
                warnings += errors
                return None

            mane = await self.mane_transcript.get_mane_transcript(
                ac, pos0, "g", end_pos=pos1, gene=classification.gene.token,
                try_longest_compatible=True
            )

            if mane:
                ac = mane["refseq"]
                pos0 = mane["pos"][0] + mane["coding_start_site"]
                pos1 = mane["pos"][1] + mane["coding_start_site"]
                classification.molecule_context = MoleculeContext.TRANSCRIPT
            else:
                return None

        outer_coords = (pos0, pos1 if pos1 else pos0)
        ival = models.SequenceInterval(
            start=models.Number(value=pos0 - 1, type="Number"),
            end=models.Number(value=pos1 if pos1 else pos0, type="Number")
        ).as_dict()

        seq_id = self.translate_sequence_identifier(ac, warnings)
        if not seq_id:
            return None

        seq_loc = self.vrs.get_sequence_loc(seq_id, ival).as_dict()

        if endpoint_name == Endpoint.NORMALIZE:
            vrs_variation = self.hgvs_dup_del_mode.interpret_variation(
                AltType.DUPLICATION, seq_loc, warnings, hgvs_dup_del_mode, ac,
                baseline_copies=baseline_copies, copy_change=copy_change
            )
        elif endpoint_name == Endpoint.HGVS_TO_COPY_NUMBER_COUNT:
            vrs_variation = self.hgvs_dup_del_mode.copy_number_count_mode(
                "dup", seq_loc, baseline_copies
            )
        elif endpoint_name == Endpoint.HGVS_TO_COPY_NUMBER_CHANGE:
            vrs_variation = self.hgvs_dup_del_mode.copy_number_change_mode(
                "dup", seq_loc, copy_change
            )
        else:
            vrs_variation = self.hgvs_dup_del_mode.default_mode(
                AltType.DUPLICATION, outer_coords, "dup", seq_loc, ac,
                baseline_copies=baseline_copies, copy_change=copy_change
            )

        if vrs_variation:
            return TranslationResult(
                vrs_variation=vrs_variation, vrs_seq_loc_ac=ac
            )
        else:
            return None

    def is_valid(
        self, gene_token, alt_ac, pos0, pos1, errors, pos2=None, pos3=None,
        residue_mode: ResidueMode = ResidueMode.RESIDUE
    ):
        """Assumes grch38"""
        gene_start = None
        gene_end = None

        for ext in gene_token.gene_descriptor.extensions:
            if ext.name == "ensembl_locations":
                if ext.value:
                    ensembl_loc = ext.value[0]
                    gene_start = ensembl_loc["interval"]["start"]["value"]
                    gene_end = ensembl_loc["interval"]["end"]["value"] - 1

        if gene_start is None and gene_end is None:
            errors.append(
                f"gene-normalizer unable to find Ensembl location for: {gene_token.token}"  # noqa: E501
            )

        for pos in [pos0, pos1, pos2, pos3]:
            if pos not in {"?", None}:
                if residue_mode == ResidueMode.RESIDUE:
                    pos -= 1

                if not (gene_start <= pos <= gene_end):
                    errors.append(
                        f"Inter-residue position {pos} out of index on {alt_ac} on gene, {gene_token.token}"  # noqa: E501
                    )
