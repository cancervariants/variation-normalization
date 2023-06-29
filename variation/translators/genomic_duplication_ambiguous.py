"""Module for Genomic Duplication Ambiguous Translation."""
from typing import Optional, List

from ga4gh.vrs import models
from ga4gh.vrsatile.pydantic.vrs_models import CopyChange

from variation.schemas.app_schemas import Endpoint
from variation.schemas.service_schema import ClinVarAssembly
from variation.schemas.token_response_schema import AltType
from variation.schemas.validation_response_schema import ValidationResult
from variation.schemas.normalize_response_schema import (
    HGVSDupDelMode as HGVSDupDelModeEnum
)
from variation.translators.translator import Translator
from variation.schemas.classification_response_schema import (
    AmbiguousType, ClassificationType, GenomicDuplicationAmbiguousClassification
)
from variation.schemas.translation_response_schema import TranslationResult
from variation.utils import get_assembly


class GenomicDuplicationAmbiguous(Translator):
    """The Genomic Duplication Ambiguous Translator class."""

    def can_translate(self, type: ClassificationType) -> bool:
        """Return if classification type is Genomic Duplication Ambiguous."""
        return type == ClassificationType.GENOMIC_DUPLICATION_AMBIGUOUS

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
        classification: GenomicDuplicationAmbiguousClassification = validation_result.classification  # noqa: E501
        vrs_variation = None
        outer_coords = None

        # TODO: Clean this up
        if do_liftover or endpoint_name == Endpoint.NORMALIZE:
            errors = []

            # Check if we need to do liftover
            assembly, w = get_assembly(self.seqrepo_access, validation_result.accession)
            if w:
                warnings.append(w)
                return None
            else:
                # assembly is either 37 or 38
                if assembly == ClinVarAssembly.GRCH37:
                    grch38_data = await self.get_grch38_data_ambiguous(
                        classification, errors, validation_result.accession
                    )
                    if errors:
                        warnings += errors
                        return None

                    ac = grch38_data["ac"]
                    pos0 = grch38_data["pos0"]
                    pos1 = grch38_data["pos1"]
                    pos2 = grch38_data["pos2"]
                    pos3 = grch38_data["pos3"]
                else:
                    ac = validation_result.accession
                    pos0 = classification.pos0
                    pos1 = classification.pos1
                    pos2 = classification.pos2
                    pos3 = classification.pos3

                assembly = ClinVarAssembly.GRCH38
        else:
            ac = validation_result.accession
            pos0 = classification.pos0
            pos1 = classification.pos1
            pos2 = classification.pos2
            pos3 = classification.pos3
            assembly = None

        if classification.gene:
            # Free text
            errors = []
            if not assembly:
                grch38_data = await self.get_grch38_data_ambiguous(
                    classification, errors, ac
                )

                if errors:
                    warnings += errors
                    return None

                self.is_valid(
                    classification.gene, grch38_data["ac"], grch38_data["pos0"],
                    grch38_data["pos1"], errors, pos2=grch38_data["pos2"],
                    pos3=grch38_data["pos3"]
                )
            else:
                self.is_valid(
                    classification.gene, ac, pos0, pos1, errors, pos2=pos2, pos3=pos3
                )

        if errors:
            warnings += errors
            return None

        if classification.ambiguous_type == AmbiguousType.AMBIGUOUS_1:
            ival = self.vrs.get_ival_certain_range(pos0, pos1, pos2, pos3)
            outer_coords = (pos0, pos3)
        elif classification.ambiguous_type == AmbiguousType.AMBIGUOUS_2:
            ival = models.SequenceInterval(
                start=self.vrs.get_start_indef_range(pos1),
                end=self.vrs.get_end_indef_range(pos2),
                type="SequenceInterval"
            ).as_dict()
            outer_coords = (pos1, pos2)
        elif classification.ambiguous_type == AmbiguousType.AMBIGUOUS_5:
            ival = models.SequenceInterval(
                start=self.vrs.get_start_indef_range(pos1),
                end=models.Number(value=pos2, type="Number"),
                type="SequenceInterval"
            ).as_dict()
            outer_coords = (pos1, pos2)
        elif classification.ambiguous_type == AmbiguousType.AMBIGUOUS_7:
            ival = models.SequenceInterval(
                start=models.Number(value=pos0 - 1, type="Number"),
                end=self.vrs.get_end_indef_range(pos2),
                type="SequenceInterval"
            ).as_dict()
            outer_coords = (pos0, pos2)

        seq_id = self.translate_sequence_identifier(ac, errors)
        if not seq_id:
            return None

        seq_loc = self.vrs.get_sequence_loc(seq_id, ival).as_dict()

        if endpoint_name == Endpoint.NORMALIZE:
            vrs_variation = self.hgvs_dup_del_mode.interpret_variation(
                AltType.DUPLICATION_AMBIGUOUS, seq_loc, warnings, hgvs_dup_del_mode, ac,
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
                AltType.DUPLICATION_AMBIGUOUS, outer_coords, "dup", seq_loc, ac,
                baseline_copies=baseline_copies, copy_change=copy_change
            )

        if vrs_variation:
            return TranslationResult(vrs_variation=vrs_variation, vrs_seq_loc_ac=ac)
        else:
            return None
