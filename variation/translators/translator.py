"""Module for translation."""
from abc import ABC, abstractmethod
from typing import List, Optional, Union

from cool_seq_tool.data_sources import MANETranscript, SeqRepoAccess, UTADatabase
from cool_seq_tool.schemas import AnnotationLayer, ResidueMode
from ga4gh.vrsatile.pydantic.vrs_models import CopyChange

from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.schemas.app_schemas import Endpoint
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.normalize_response_schema import HGVSDupDelModeOption
from variation.schemas.token_response_schema import AltType, GeneToken
from variation.schemas.translation_response_schema import (
    TranslationResult,
    VrsSeqLocAcStatus,
)
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.genomic_base import GenomicBase
from variation.vrs_representation import VRSRepresentation


class Translator(ABC):
    """Class for translating to VRS representations"""

    def __init__(
        self,
        seqrepo_access: SeqRepoAccess,
        mane_transcript: MANETranscript,
        uta: UTADatabase,
        vrs: VRSRepresentation,
        hgvs_dup_del_mode: HGVSDupDelMode,
    ) -> None:
        """Initialize the Translator class.

        :param seqrepo_access: Access to SeqRepo data
        :param mane_transcript: Access MANE Transcript information
        :param uta: Access to UTA queries
        :param vrs: Class for creating VRS objects
        :param hgvs_dup_del_mode: Class for interpreting HGVS duplications and deletions
        """
        self.seqrepo_access = seqrepo_access
        self.uta = uta
        self.genomic_base = GenomicBase(self.seqrepo_access, self.uta)
        self.mane_transcript = mane_transcript
        self.vrs = vrs
        self.hgvs_dup_del_mode = hgvs_dup_del_mode

    @abstractmethod
    def can_translate(self, classification_type: ClassificationType) -> bool:
        """Determine if it's possible to translate a classification.

        :param classification_type: Classification type found
        :return: `True` if `classification_type` matches translator's classification
            type. Otherwise, `False`
        """

    @abstractmethod
    async def translate(
        self,
        validation_result: ValidationResult,
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

    def translate_sequence_identifier(
        self, sequence_id: str, errors: List[str]
    ) -> Optional[str]:
        """Translate `sequence_id` to ga4gh identifier

        :param sequence_id: Sequence ID to translate
        :param errors: List of errors. This will get mutated if an error occurs when
            attempting to get ga4gh identifier
        :return: GA4GH Sequence Identifier if successful, else `None`
        """
        ga4gh_seq_id = None
        try:
            ids = self.seqrepo_access.translate_sequence_identifier(
                sequence_id, "ga4gh"
            )
        except KeyError as e:
            errors.append(str(e))
        else:
            if not ids:
                errors.append(
                    f"Unable to find ga4gh sequence identifiers for: {sequence_id}"
                )

            ga4gh_seq_id = ids[0]
        return ga4gh_seq_id

    def is_valid(
        self,
        gene_token: GeneToken,
        alt_ac: str,
        pos0: int,
        pos1: int,
        errors: List[str],
        pos2: Optional[int] = None,
        pos3: Optional[int] = None,
        residue_mode: ResidueMode = ResidueMode.RESIDUE,
    ) -> None:
        """Check that positions are valid on a gene. Will mutate `errors` if invalid.

        :param gene_token: Gene token
        :param alt_ac: Genomic RefSeq accession
        :param pos0: Position 0 (GRCh38 assembly)
        :param pos1: Position 1 (GRCh38 assembly)
        :param errors: List of errors. Will get mutated if invalid.
        :param pos2: Position 2 (GRCh38 assembly)
        :param pos3: Position 3 (GRCh38 assembly)
        :param residue_mode: Residue mode for positions
        """
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

    async def get_p_or_cdna_translation_result(
        self,
        endpoint_name: Endpoint,
        validation_result: ValidationResult,
        start_pos: int,
        end_pos: int,
        alt_type: AltType,
        coordinate_type: Union[AnnotationLayer.PROTEIN, AnnotationLayer.CDNA],
        errors: List[str],
        cds_start: Optional[int] = None,
        ref: Optional[str] = None,
        alt: Optional[str] = None,
    ) -> Optional[TranslationResult]:
        """Get translation result for validation result. Used for unambiguous
        variations on protein or cDNA coordinate types

        :param endpoint_name: Name of endpoint that is being used
        :param validation_result: Validation result for a classification
        :param start_pos: Start position (residue-mode)
        :param end_pos: End position (residue-mode)
        :param alt_type: Alteration type for validation result
        :param coordinate_type: Coordinate type for validation result
        :param errors: List of errors. Will be mutated if errors are found
        :param cds_start: Coding start site. Only required for
            `coordinate_type == AnnotationLayer.CDNA`.
        :param ref: Expected reference sequence
        :param alt: Expected change
        :return: Translation result if successful. Else, `None`
        """
        vrs_allele = None
        vrs_seq_loc_ac = None
        vrs_seq_loc_ac_status = VrsSeqLocAcStatus.NA

        if endpoint_name == Endpoint.NORMALIZE:
            mane = await self.mane_transcript.get_mane_transcript(
                validation_result.accession,
                start_pos,
                coordinate_type,
                end_pos=end_pos,
                try_longest_compatible=True,
                residue_mode=ResidueMode.RESIDUE.value,
                ref=ref,
            )

            if mane:
                vrs_seq_loc_ac = mane["refseq"]
                vrs_seq_loc_ac_status = mane["status"]
                vrs_allele = self.vrs.to_vrs_allele(
                    vrs_seq_loc_ac,
                    mane["pos"][0] + 1,
                    mane["pos"][1] + 1,
                    coordinate_type,
                    alt_type,
                    errors,
                    cds_start=mane.get("coding_start_site", None),
                    alt=alt,
                )

        if not vrs_allele:
            vrs_seq_loc_ac = validation_result.accession
            vrs_allele = self.vrs.to_vrs_allele(
                vrs_seq_loc_ac,
                start_pos,
                end_pos,
                coordinate_type,
                alt_type,
                errors,
                cds_start=cds_start,
                alt=alt,
            )

        if vrs_allele and vrs_seq_loc_ac:
            return TranslationResult(
                vrs_variation=vrs_allele,
                vrs_seq_loc_ac=vrs_seq_loc_ac,
                vrs_seq_loc_ac_status=vrs_seq_loc_ac_status,
                og_ac=validation_result.accession,
                validation_result=validation_result,
            )
        else:
            return None
