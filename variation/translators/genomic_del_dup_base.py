"""Module for Genomic Deletion Translation."""
from typing import List, NamedTuple, Optional, Union

from cool_seq_tool.schemas import ResidueMode
from ga4gh.vrs import models
from ga4gh.vrsatile.pydantic.vrs_models import CopyChange
from ga4gh.vrsatile.pydantic.vrsatile_models import MoleculeContext
from pydantic import StrictInt, StrictStr, ValidationError

from variation.schemas.app_schemas import Endpoint
from variation.schemas.classification_response_schema import (
    GenomicDeletionClassification,
    GenomicDuplicationClassification,
    Nomenclature,
)
from variation.schemas.normalize_response_schema import HGVSDupDelModeOption
from variation.schemas.service_schema import ClinVarAssembly
from variation.schemas.token_response_schema import AltType
from variation.schemas.translation_response_schema import (
    TranslationResult,
    VrsSeqLocAcStatus,
)
from variation.schemas.validation_response_schema import ValidationResult
from variation.translators.translator import Translator
from variation.utils import get_assembly


class DelDupData(NamedTuple):
    """Represents genomic dup/del data"""

    ac: StrictStr
    pos0: StrictInt
    pos1: Optional[StrictInt]


class GenomicDelDupTranslator(Translator):
    """Class for translating genomic deletions and duplications to VRS
    representations
    """

    async def get_grch38_data(
        self,
        classification: Union[
            GenomicDeletionClassification, GenomicDuplicationClassification
        ],
        errors: List[str],
        ac: str,
    ) -> DelDupData:
        """Get GRCh38 data for genomic duplication or deletion classification

        :param classification: Classification to get translation for
        :param errors: List of errors. Will be mutated if errors are found
        :param ac: Genomic RefSeq accession
        :return: Data on GRCh38 assembly if successful liftover. Else, `None`
        """
        pos0, pos1, new_ac = None, None, None

        if classification.pos1:
            grch38_pos = await self.mane_transcript.g_to_grch38(
                ac, classification.pos0, classification.pos1
            )
            if grch38_pos:
                pos0, pos1 = grch38_pos["pos"]
                new_ac = grch38_pos["ac"]
        else:
            grch38_pos = await self.mane_transcript.g_to_grch38(
                ac, classification.pos0, classification.pos0
            )
            if grch38_pos:
                pos0, _ = grch38_pos["pos"]
                new_ac = grch38_pos["ac"]

        if not new_ac:
            errors.append(f"Unable to find a GRCh38 accession for: {ac}")

        try:
            data = DelDupData(ac=new_ac, pos0=pos0, pos1=pos1)
        except ValidationError:
            data = None
        return data

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
        if isinstance(classification, GenomicDeletionClassification):
            alt_type = AltType.DELETION
        else:
            alt_type = AltType.DUPLICATION

        grch38_data = None
        vrs_variation = None
        vrs_seq_loc_ac_status = VrsSeqLocAcStatus.NA

        if do_liftover or endpoint_name == Endpoint.NORMALIZE:
            errors = []
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

                    pos0 = grch38_data.pos0
                    pos1 = grch38_data.pos1
                    ac = grch38_data.ac
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

        if all(
            (
                endpoint_name == Endpoint.NORMALIZE,
                classification.nomenclature == Nomenclature.FREE_TEXT,
                classification.gene_token,
            )
        ):
            errors = []
            if not assembly and not grch38_data:
                grch38_data = await self.get_grch38_data(classification, errors, ac)
                if errors:
                    warnings += errors
                    return None

                pos0 = grch38_data.pos0
                pos1 = grch38_data.pos1
                ac = grch38_data.ac

            self.is_valid(classification.gene_token, ac, pos0, pos1, errors)

            if errors:
                warnings += errors
                return None

            mane = await self.mane_transcript.get_mane_transcript(
                ac,
                pos0,
                "g",
                end_pos=pos1,
                try_longest_compatible=True,
                residue_mode=ResidueMode.RESIDUE,
                gene=classification.gene_token.token
                if classification.gene_token
                else None,
            )

            if mane:
                # mane is 0 - based, but we are using residue
                ac = mane["refseq"]
                vrs_seq_loc_ac_status = mane["status"]
                pos0 = mane["pos"][0] + mane["coding_start_site"] + 1
                pos1 = mane["pos"][1] + mane["coding_start_site"] + 1
                classification.molecule_context = MoleculeContext.TRANSCRIPT
            else:
                return None

        outer_coords = (pos0, pos1 if pos1 else pos0)

        ival = models.SequenceInterval(
            start=models.Number(value=pos0 - 1, type="Number"),
            end=models.Number(value=pos1 if pos1 else pos0, type="Number"),
        ).as_dict()

        seq_id = self.translate_sequence_identifier(ac, warnings)
        if not seq_id:
            return None

        seq_loc = self.vrs.get_sequence_loc(seq_id, ival).as_dict()

        if endpoint_name == Endpoint.NORMALIZE:
            vrs_variation = self.hgvs_dup_del_mode.interpret_variation(
                alt_type,
                seq_loc,
                warnings,
                hgvs_dup_del_mode,
                ac,
                pos=outer_coords,
                baseline_copies=baseline_copies,
                copy_change=copy_change,
            )
        elif endpoint_name == Endpoint.HGVS_TO_COPY_NUMBER_COUNT:
            vrs_variation = self.hgvs_dup_del_mode.copy_number_count_mode(
                alt_type, seq_loc, baseline_copies
            )
        elif endpoint_name == Endpoint.HGVS_TO_COPY_NUMBER_CHANGE:
            vrs_variation = self.hgvs_dup_del_mode.copy_number_change_mode(
                alt_type, seq_loc, copy_change
            )
        else:
            vrs_variation = self.hgvs_dup_del_mode.default_mode(
                alt_type,
                outer_coords,
                seq_loc,
                ac,
                baseline_copies=baseline_copies,
                copy_change=copy_change,
            )

        if vrs_variation:
            return TranslationResult(
                vrs_variation=vrs_variation,
                vrs_seq_loc_ac=ac,
                vrs_seq_loc_ac_status=vrs_seq_loc_ac_status,
                og_ac=validation_result.accession,
                validation_result=validation_result,
            )
        else:
            return None
