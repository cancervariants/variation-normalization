"""Module for translating genomic ambiguous deletions and duplications"""
from typing import Dict, List, Literal, NamedTuple, Optional, Tuple, Union

from ga4gh.vrs import models
from ga4gh.vrsatile.pydantic.vrs_models import CopyChange
from pydantic import StrictInt, StrictStr, ValidationError

from variation.schemas.app_schemas import Endpoint
from variation.schemas.classification_response_schema import (
    AmbiguousType,
    GenomicDeletionAmbiguousClassification,
    GenomicDuplicationAmbiguousClassification,
    Nomenclature,
)
from variation.schemas.normalize_response_schema import HGVSDupDelModeOption
from variation.schemas.service_schema import ClinVarAssembly
from variation.schemas.token_response_schema import AltType
from variation.schemas.translation_response_schema import TranslationResult
from variation.schemas.validation_response_schema import ValidationResult
from variation.translators.translator import Translator
from variation.utils import get_assembly


class AmbiguousData(NamedTuple):
    """Represents Ambiguous data"""

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
    """Class for translating genomic ambiguous deletions and duplications to VRS
    representations
    """

    async def get_grch38_data_ambiguous(
        self,
        classification: Union[
            GenomicDeletionAmbiguousClassification,
            GenomicDuplicationAmbiguousClassification,
        ],
        errors: List[str],
        ac: str,
    ) -> Optional[AmbiguousData]:
        """Get GRCh38 data for genomic ambiguous duplication or deletion classification

        :param classification: Classification to get translation for
        :param errors: List of errors. Will be mutated if errors are found
        :param ac: Genomic RefSeq accession
        :return: Ambiguous data on GRCh38 assembly if successful liftover. Else, `None`
        """
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
        elif classification.ambiguous_type in {
            AmbiguousType.AMBIGUOUS_2,
            AmbiguousType.AMBIGUOUS_5,
        }:
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
                (
                    pos0,
                    pos2,
                ) = grch38["pos"]
                new_ac = grch38["ac"]

        if not new_ac:
            errors.append(f"Unable to find a GRCh38 accession for: {ac}")

        try:
            ambiguous_data = AmbiguousData(
                ac=new_ac, pos0=pos0, pos1=pos1, pos2=pos2, pos3=pos3
            )
        except ValidationError:
            ambiguous_data = None
        return ambiguous_data

    def get_dup_del_ambiguous_seq_loc(
        self,
        ambiguous_type: AmbiguousType,
        ac: str,
        pos0: Union[int, Literal["?"]],
        pos1: Optional[Union[int, Literal["?"]]],
        pos2: Union[int, Literal["?"]],
        pos3: Optional[Union[int, Literal["?"]]],
        warnings: List[str],
    ) -> AmbiguousSequenceLocation:
        """Get VRS Sequence Location and outer coords for genomic ambiguous duplication
        or deletion

        :param ambiguous_type: Type of ambiguous expression used
        :param ac: Genomic RefSeq accession
        :param pos0: Position 0 (residue)
        :param pos1: Position 1 (residue)
        :param pos2: Position 2 (residue)
        :param pos3: Position 3 (residue)
        :param warnings: List of warnings
        :return: Ambiguous sequence location data containing VRS Sequence Location
            and outer coordinates
        """
        if ambiguous_type == AmbiguousType.AMBIGUOUS_1:
            ival = self.vrs.get_ival_certain_range(pos0, pos1, pos2, pos3)
            outer_coords = (pos0, pos3)
        elif ambiguous_type == AmbiguousType.AMBIGUOUS_2:
            ival = models.SequenceInterval(
                start=self.vrs.get_start_indef_range(pos1),
                end=self.vrs.get_end_indef_range(pos2),
                type="SequenceInterval",
            ).as_dict()
            outer_coords = (pos1, pos2)
        elif ambiguous_type == AmbiguousType.AMBIGUOUS_5:
            ival = models.SequenceInterval(
                start=self.vrs.get_start_indef_range(pos1),
                end=models.Number(value=pos2, type="Number"),
                type="SequenceInterval",
            ).as_dict()
            outer_coords = (pos1, pos2)
        elif ambiguous_type == AmbiguousType.AMBIGUOUS_7:
            ival = models.SequenceInterval(
                start=models.Number(value=pos0 - 1, type="Number"),
                end=self.vrs.get_end_indef_range(pos2),
                type="SequenceInterval",
            ).as_dict()
            outer_coords = (pos0, pos2)
        # No else since validator should catch if the ambiguous type is supported or not

        seq_id = self.translate_sequence_identifier(ac, warnings)

        return AmbiguousSequenceLocation(
            seq_loc=self.vrs.get_sequence_loc(seq_id, ival).as_dict()
            if seq_id
            else None,
            outer_coords=outer_coords,
        )

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
        # First will translate valid result to VRS Allele
        classification = validation_result.classification
        if isinstance(classification, GenomicDeletionAmbiguousClassification):
            alt_type = AltType.DELETION_AMBIGUOUS
        else:
            alt_type = AltType.DUPLICATION_AMBIGUOUS

        grch38_data = None

        if do_liftover or endpoint_name == Endpoint.NORMALIZE:
            errors = []
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

        if all(
            (
                endpoint_name == Endpoint.NORMALIZE,
                classification.nomenclature == Nomenclature.FREE_TEXT,
                classification.gene_token,
            )
        ):
            errors = []
            if not assembly and not grch38_data:
                grch38_data = await self.get_grch38_data_ambiguous(
                    classification, errors, ac
                )

                if errors:
                    warnings += errors
                    return None

                ac = grch38_data.ac
                pos0 = grch38_data.pos0
                pos1 = grch38_data.pos1
                pos2 = grch38_data.pos2
                pos3 = grch38_data.pos3

                self.is_valid(
                    classification.gene_token,
                    ac,
                    pos0,
                    pos1,
                    errors,
                    pos2=pos2,
                    pos3=pos3,
                )
            else:
                self.is_valid(
                    classification.gene_token,
                    ac,
                    pos0,
                    pos1,
                    errors,
                    pos2=pos2,
                    pos3=pos3,
                )

            if errors:
                warnings += errors
                return None

        ambiguous_seq_loc_data = self.get_dup_del_ambiguous_seq_loc(
            classification.ambiguous_type, ac, pos0, pos1, pos2, pos3, warnings
        )

        if endpoint_name == Endpoint.NORMALIZE:
            vrs_variation = self.hgvs_dup_del_mode.interpret_variation(
                alt_type,
                ambiguous_seq_loc_data.seq_loc,
                warnings,
                hgvs_dup_del_mode,
                ac,
                pos=ambiguous_seq_loc_data.outer_coords,
                baseline_copies=baseline_copies,
                copy_change=copy_change,
            )
        elif endpoint_name == Endpoint.HGVS_TO_COPY_NUMBER_COUNT:
            vrs_variation = self.hgvs_dup_del_mode.copy_number_count_mode(
                alt_type, ambiguous_seq_loc_data.seq_loc, baseline_copies
            )
        elif endpoint_name == Endpoint.HGVS_TO_COPY_NUMBER_CHANGE:
            vrs_variation = self.hgvs_dup_del_mode.copy_number_change_mode(
                alt_type, ambiguous_seq_loc_data.seq_loc, copy_change
            )
        else:
            vrs_variation = self.hgvs_dup_del_mode.default_mode(
                alt_type,
                ambiguous_seq_loc_data.outer_coords,
                ambiguous_seq_loc_data.seq_loc,
                ac,
                baseline_copies=baseline_copies,
                copy_change=copy_change,
            )

        if vrs_variation:
            return TranslationResult(
                vrs_variation=vrs_variation,
                vrs_seq_loc_ac=ac,
                og_ac=validation_result.accession,
                validation_result=validation_result,
            )
        else:
            return None
