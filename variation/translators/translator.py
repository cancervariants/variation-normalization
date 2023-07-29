"""Module for translation."""
from abc import ABC, abstractmethod
from typing import List, Optional

from cool_seq_tool.data_sources import MANETranscript, SeqRepoAccess, UTADatabase
from cool_seq_tool.schemas import ResidueMode
from ga4gh.vrsatile.pydantic.vrs_models import CopyChange

from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.schemas.app_schemas import Endpoint
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.normalize_response_schema import HGVSDupDelModeOption
from variation.schemas.token_response_schema import GeneToken
from variation.schemas.translation_response_schema import TranslationResult
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

    def validate_reference_sequence(
        self,
        ac: str,
        start_pos: int,
        end_pos: int,
        expected_ref: str,
        residue_mode: ResidueMode = ResidueMode.RESIDUE,
    ) -> Optional[str]:
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
            ac, start=start_pos, end=end_pos, residue_mode=residue_mode
        )

        if not err_msg and (actual_ref != expected_ref):
            err_msg = (
                f"Expected to find {expected_ref} at positions ({start_pos}, "
                f"{end_pos}) on {ac} but found {actual_ref}"
            )

        return err_msg
