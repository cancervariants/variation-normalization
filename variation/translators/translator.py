"""Module for translation."""
from abc import abstractmethod, ABC
from typing import Optional, List

from ga4gh.vrsatile.pydantic.vrs_models import CopyChange
from gene.query import QueryHandler as GeneQueryHandler
from cool_seq_tool.data_sources import (
    TranscriptMappings, SeqRepoAccess, UTADatabase, MANETranscript
)
from cool_seq_tool.schemas import ResidueMode

from variation.schemas.token_response_schema import GeneToken
from variation.validators.genomic_base import GenomicBase
from variation.vrs_representation import VRSRepresentation
from variation.schemas.app_schemas import Endpoint
from variation.schemas.validation_response_schema import ValidationResult
from variation.schemas.normalize_response_schema import (
    HGVSDupDelMode as HGVSDupDelModeEnum
)
from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.schemas.classification_response_schema import ClassificationType
from variation.schemas.translation_response_schema import TranslationResult


class Translator(ABC):
    """The translation class."""

    def __init__(
        self,
        seqrepo_access: SeqRepoAccess,
        mane_transcript: MANETranscript,
        uta: UTADatabase,
        vrs: VRSRepresentation,
        hgvs_dup_del_mode: HGVSDupDelMode
    ) -> None:
        """Initialize the DelIns validator.

        :param seqrepo_access: Access to SeqRepo data
        :param gene_symbol: Gene symbol tokenizer
        :param mane_transcript: Access MANE Transcript information
        :param uta: Access to UTA queries
        :param vrs: Class for creating VRS objects
        """
        self.seqrepo_access = seqrepo_access
        self.uta = uta
        self.genomic_base = GenomicBase(self.seqrepo_access, self.uta)
        self.mane_transcript = mane_transcript
        self.vrs = vrs
        self.hgvs_dup_del_mode = hgvs_dup_del_mode

    @abstractmethod
    def can_translate(self, type: ClassificationType) -> bool:
        """Determine if it's possible to translate a classification."""

    @abstractmethod
    async def translate(
        self,
        validation_result: ValidationResult,
        endpoint_name: Optional[Endpoint] = None,
        hgvs_dup_del_mode: HGVSDupDelModeEnum = HGVSDupDelModeEnum.DEFAULT,
        baseline_copies: Optional[int] = None,
        copy_change: Optional[CopyChange] = None,
        do_liftover: bool = False
    ) -> Optional[TranslationResult]:
        """"""

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
        self, gene_token: GeneToken, alt_ac: str, pos0: int, pos1: int,
        errors: List[str], pos2: Optional[int] = None, pos3: Optional[int] = None,
        residue_mode: ResidueMode = ResidueMode.RESIDUE
    ) -> None:
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

    def validate_reference_sequence(
        self, ac: str, start_pos: int, end_pos: int,
        expected_ref: str, residue_mode: ResidueMode = ResidueMode.RESIDUE
    ) -> Optional[str]:
        """Validate that expected reference sequence matches actual"""
        actual_ref, _ = self.seqrepo_access.get_reference_sequence(
            ac, start=start_pos, end=end_pos, residue_mode=residue_mode
        )

        msg = None
        if actual_ref != expected_ref:
            msg = (f"Expected to find {expected_ref} at positions "
                   f"({start_pos}, {end_pos}) on {ac} but found {actual_ref}")

        return msg
