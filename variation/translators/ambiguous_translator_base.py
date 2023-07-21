from typing import Dict, Optional, List, NamedTuple, Tuple, Union, Literal

from ga4gh.vrs import models
from ga4gh.vrsatile.pydantic.vrs_models import CopyChange
from pydantic import StrictInt, StrictStr

from variation.schemas.app_schemas import Endpoint
from variation.schemas.validation_response_schema import ValidationResult
from variation.schemas.normalize_response_schema import (
    HGVSDupDelMode as HGVSDupDelModeEnum
)
from variation.schemas.classification_response_schema import (
    AmbiguousType, GenomicDeletionAmbiguousClassification
)
from variation.schemas.token_response_schema import AltType
from variation.schemas.translation_response_schema import TranslationResult
from variation.translators.translator import Translator
from variation.schemas.service_schema import ClinVarAssembly
from variation.utils import get_assembly


class AmbiguousData(NamedTuple):
    """Represents Ambiguous data assembly"""

    ac: StrictStr
    pos0: Union[StrictInt, Literal["?"]]
    pos1: Optional[Union[StrictInt, Literal["?"]]]
    pos2: Union[StrictInt, Literal["?"]]
    pos3: Optional[Union[StrictInt, Literal["?"]]]


class AmbiguousSequenceLocation(NamedTuple):
    """Represents ambiguous dup/del sequence location info"""

    seq_loc: Optional[Dict]
    outer_coords: Tuple[int, int]


class AmbiguousTranslator(Translator):

    async def get_grch38_data_ambiguous(
        self, classification, errors, ac
    ) -> AmbiguousData:
        pos0, pos1, pos2, pos3, new_ac = None, None, None, None, None
        if classification.ambiguous_type == AmbiguousType.AMBIGUOUS_1:
            grch38_pos0_pos1 = await self.mane_transcript.g_to_grch38(
                ac, classification.pos0, classification.pos1
            )
            if grch38_pos0_pos1:
                pos0, pos1 = grch38_pos0_pos1["pos"]
                ac_pos0_pos1 = grch38_pos0_pos1["ac"]

                grch38_pos2_pos3 = await self.mane_transcript.g_to_grch38(
                    ac, classification.pos2, classification.pos3
                )

                if grch38_pos2_pos3:
                    pos2, pos3 = grch38_pos2_pos3["pos"]
                    ac_pos2_pos3 = grch38_pos2_pos3["ac"]

                    if ac_pos0_pos1 != ac_pos2_pos3:
                        errors.append(
                            f"{ac_pos0_pos1} does not equal {ac_pos2_pos3} when lifting"
                            " over to GRCh38"
                        )
                    else:
                        new_ac = ac_pos0_pos1
        elif classification.ambiguous_type in {AmbiguousType.AMBIGUOUS_2,
                                               AmbiguousType.AMBIGUOUS_5}:
            grch38 = await self.mane_transcript.g_to_grch38(
                ac, classification.pos1, classification.pos2
            )
            if grch38:
                pos1, pos2 = grch38["pos"]
                new_ac = grch38["ac"]
        elif classification.ambiguous_type == AmbiguousType.AMBIGUOUS_7:
            grch38 = await self.mane_transcript.g_to_grch38(
                ac, classification.pos0, classification.pos2
            )
            if grch38:
                pos0, pos2, = grch38["pos"]
                new_ac = grch38["ac"]

        if not new_ac:
            errors.append(f"Unable to find a GRCh38 accession for: {ac}")

        return AmbiguousData(ac=new_ac, pos0=pos0, pos1=pos1, pos2=pos2, pos3=pos3)

    def get_dup_del_ambiguous_seq_loc(
        self, ambiguous_type: AmbiguousType, ac: str, pos0, pos1, pos2, pos3,
        warnings: List[str]
    ) -> AmbiguousSequenceLocation:
        if ambiguous_type == AmbiguousType.AMBIGUOUS_1:
            ival = self.vrs.get_ival_certain_range(pos0, pos1, pos2, pos3)
            outer_coords = (pos0, pos3)
        elif ambiguous_type == AmbiguousType.AMBIGUOUS_2:
            ival = models.SequenceInterval(
                start=self.vrs.get_start_indef_range(pos1),
                end=self.vrs.get_end_indef_range(pos2),
                type="SequenceInterval"
            ).as_dict()
            outer_coords = (pos1, pos2)
        elif ambiguous_type == AmbiguousType.AMBIGUOUS_5:
            ival = models.SequenceInterval(
                start=self.vrs.get_start_indef_range(pos1),
                end=models.Number(value=pos2, type="Number"),
                type="SequenceInterval"
            ).as_dict()
            outer_coords = (pos1, pos2)
        elif ambiguous_type == AmbiguousType.AMBIGUOUS_7:
            ival = models.SequenceInterval(
                start=models.Number(value=pos0 - 1, type="Number"),
                end=self.vrs.get_end_indef_range(pos2),
                type="SequenceInterval"
            ).as_dict()
            outer_coords = (pos0, pos2)

        seq_id = self.translate_sequence_identifier(ac, warnings)

        return AmbiguousSequenceLocation(
            seq_loc=self.vrs.get_sequence_loc(seq_id, ival).as_dict() if seq_id else None,  # noqa: E501
            outer_coords=outer_coords
        )

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
        classification = validation_result.classification
        if isinstance(classification, GenomicDeletionAmbiguousClassification):
            alt_type = AltType.DELETION_AMBIGUOUS
            del_or_dup = "del"
        else:
            alt_type = AltType.DUPLICATION_AMBIGUOUS
            del_or_dup = "dup"

        grch38_data = None

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

                    ac = grch38_data.ac
                    pos0 = grch38_data.pos0
                    pos1 = grch38_data.pos1
                    pos2 = grch38_data.pos2
                    pos3 = grch38_data.pos3
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

        if classification.gene_token:
            # Free text
            errors = []
            if not assembly:
                if not grch38_data:
                    grch38_data = await self.get_grch38_data_ambiguous(
                        classification, errors, ac
                    )

                    if errors:
                        warnings += errors
                        return None

                self.is_valid(
                    classification.gene_token, grch38_data.ac, grch38_data.pos0,
                    grch38_data.pos1, errors, pos2=grch38_data.pos2,
                    pos3=grch38_data.pos3
                )
            else:
                self.is_valid(
                    classification.gene_token, ac, pos0, pos1, errors, pos2=pos2,
                    pos3=pos3
                )

            if errors:
                warnings += errors
                return None

        ambiguous_seq_loc_data = self.get_dup_del_ambiguous_seq_loc(
            classification.ambiguous_type, ac, pos0, pos1, pos2, pos3, warnings
        )

        if endpoint_name == Endpoint.NORMALIZE:
            vrs_variation = self.hgvs_dup_del_mode.interpret_variation(
                alt_type, ambiguous_seq_loc_data.seq_loc, warnings,
                hgvs_dup_del_mode, ac, baseline_copies=baseline_copies,
                copy_change=copy_change
            )
        elif endpoint_name == Endpoint.HGVS_TO_COPY_NUMBER_COUNT:
            vrs_variation = self.hgvs_dup_del_mode.copy_number_count_mode(
                del_or_dup, ambiguous_seq_loc_data.seq_loc, baseline_copies
            )
        elif endpoint_name == Endpoint.HGVS_TO_COPY_NUMBER_CHANGE:
            vrs_variation = self.hgvs_dup_del_mode.copy_number_change_mode(
                del_or_dup, ambiguous_seq_loc_data.seq_loc, copy_change
            )
        else:
            vrs_variation = self.hgvs_dup_del_mode.default_mode(
                alt_type, ambiguous_seq_loc_data.outer_coords,
                del_or_dup, ambiguous_seq_loc_data.seq_loc, ac,
                baseline_copies=baseline_copies, copy_change=copy_change
            )

        if vrs_variation:
            return TranslationResult(
                vrs_variation=vrs_variation, vrs_seq_loc_ac=ac,
                og_ac=validation_result.accession, validation_result=validation_result
            )
        else:
            return None
