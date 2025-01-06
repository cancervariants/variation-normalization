"""Module for translation."""

from abc import ABC, abstractmethod

from cool_seq_tool.handlers import SeqRepoAccess
from cool_seq_tool.mappers import ManeTranscript
from cool_seq_tool.schemas import AnnotationLayer, CoordinateType, ManeGeneData
from cool_seq_tool.sources import UtaDatabase
from ga4gh.core.models import Extension
from ga4gh.vrs import models

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
        mane_transcript: ManeTranscript,
        uta: UtaDatabase,
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
        endpoint_name: Endpoint | None = None,
        hgvs_dup_del_mode: HGVSDupDelModeOption = HGVSDupDelModeOption.DEFAULT,
        baseline_copies: int | None = None,
        copy_change: models.CopyChange | None = None,
        do_liftover: bool = False,
    ) -> TranslationResult | None:
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

    def is_valid(
        self,
        gene_token: GeneToken,
        alt_ac: str,
        pos0: int,
        pos1: int,
        errors: list[str],
        pos2: int | None = None,
        pos3: int | None = None,
        residue_mode: CoordinateType = CoordinateType.RESIDUE,
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

        for ext in gene_token.gene.extensions:
            if ext.name == "ensembl_locations" and ext.value:
                ensembl_loc = ext.value[0]
                gene_start = ensembl_loc["start"]
                gene_end = ensembl_loc["end"] - 1

        if gene_start is None and gene_end is None:
            errors.append(
                f"gene-normalizer unable to find Ensembl location for: {gene_token.token}"
            )

        for pos in [pos0, pos1, pos2, pos3]:
            if pos not in {"?", None}:
                if residue_mode == CoordinateType.RESIDUE:
                    pos -= 1

                if not (gene_start <= pos <= gene_end):
                    errors.append(
                        f"Inter-residue position {pos} out of index on {alt_ac} on gene, {gene_token.token}"
                    )

    def validate_reference_sequence(
        self,
        ac: str,
        start_pos: int,
        end_pos: int,
        expected_ref: str,
        residue_mode: CoordinateType = CoordinateType.RESIDUE,
    ) -> str | None:
        """Validate that expected reference sequence matches actual reference sequence
        This is also in validator, but there is a ticket to have this method be moved
        to cool-seq-tool. Once added, will be removed

        :param ac: Accession
        :param start_pos: Start position
        :param end_pos: End position
        :param expected_ref: The expected reference sequence (from input string)
        :param residue_mode: Residue mode for `start_pos` and `end_pos`
        :return: Invalid message if invalid. If valid, `None`
        """
        actual_ref, err_msg = self.seqrepo_access.get_reference_sequence(
            ac, start=start_pos, end=end_pos, coordinate_type=residue_mode
        )

        if not err_msg and (actual_ref != expected_ref):
            err_msg = (
                f"Expected to find {expected_ref} at positions ({start_pos}, "
                f"{end_pos}) on {ac} but found {actual_ref}"
            )

        return err_msg

    async def get_p_or_cdna_translation_result(
        self,
        endpoint_name: Endpoint,
        validation_result: ValidationResult,
        start_pos: int,
        end_pos: int,
        alt_type: AltType,
        coordinate_type: AnnotationLayer,
        errors: list[str],
        cds_start: int | None = None,
        ref: str | None = None,
        alt: str | None = None,
    ) -> TranslationResult | None:
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
        :raises ValueError: If ``coordinate`` type not one of
            ``AnnotationLayer.PROTEIN`` or ``AnnotationLayer.CDNA``
        :return: Translation result if successful. Else, `None`
        """
        supported_coordinate_types = {AnnotationLayer.PROTEIN, AnnotationLayer.CDNA}
        if coordinate_type not in supported_coordinate_types:
            err_msg = f"`coordinate_type` must be one of {supported_coordinate_types}"
            raise ValueError(err_msg)

        vrs_allele = None
        vrs_seq_loc_ac = None
        vrs_seq_loc_ac_status = VrsSeqLocAcStatus.NA

        if endpoint_name == Endpoint.NORMALIZE:
            mane = await self.mane_transcript.get_mane_transcript(
                validation_result.accession,
                start_pos,
                end_pos if end_pos is not None else start_pos,
                coordinate_type,
                try_longest_compatible=True,
                coordinate_type=CoordinateType.RESIDUE,
                ref=ref,
            )

            if mane:
                vrs_seq_loc_ac = mane.refseq
                vrs_seq_loc_ac_status = mane.status

                try:
                    cds_start = mane.coding_start_site
                except AttributeError:
                    cds_start = None

                vrs_allele = self.vrs.to_vrs_allele(
                    vrs_seq_loc_ac,
                    mane.pos[0],
                    mane.pos[1],
                    coordinate_type,
                    alt_type,
                    errors,
                    cds_start=cds_start,
                    alt=alt,
                    residue_mode=CoordinateType.INTER_RESIDUE,
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
                residue_mode=CoordinateType.RESIDUE,
            )

        if vrs_allele and vrs_seq_loc_ac:
            return TranslationResult(
                vrs_variation=vrs_allele,
                vrs_seq_loc_ac=vrs_seq_loc_ac,
                vrs_seq_loc_ac_status=vrs_seq_loc_ac_status,
                og_ac=validation_result.accession,
                validation_result=validation_result,
            )

        return None

    @staticmethod
    def _mane_gene_extensions(
        mane_genes: list[ManeGeneData],
    ) -> list[Extension] | None:
        """Transform mane genes to list of extensions

        This is only used in Genomic translators

        :param mane_genes: Optional list of mane gene data
        :return: List of extensions containing mane gene data if found. Otherwise,
            ``None``
        """
        mane_genes_exts = None
        if mane_genes:
            mane_genes_exts = [
                Extension(
                    name="mane_genes",
                    value=mane_genes,
                )
            ]
        return mane_genes_exts
