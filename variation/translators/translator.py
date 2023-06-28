"""Module for translation."""
from abc import abstractmethod, ABC
from typing import Dict, Optional, List

from ga4gh.vrsatile.pydantic.vrs_models import CopyChange
from gene.query import QueryHandler as GeneQueryHandler
from cool_seq_tool.data_sources import (
    TranscriptMappings, SeqRepoAccess, UTADatabase, MANETranscript
)

from variation.validators.genomic_base import GenomicBase
from variation.tokenizers import GeneSymbol
from variation.vrs_representation import VRSRepresentation
from variation.schemas.app_schemas import Endpoint
from variation.schemas.validation_response_schema import ValidationResult
from variation.schemas.normalize_response_schema import (
    HGVSDupDelMode as HGVSDupDelModeEnum
)
from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.schemas.classification_response_schema import AmbiguousType, ClassificationType
from variation.schemas.translation_response_schema import TranslationResult


class Translator(ABC):
    """The translation class."""

    def __init__(
        self,
        seqrepo_access: SeqRepoAccess,
        transcript_mappings: TranscriptMappings,
        gene_symbol: GeneSymbol,
        mane_transcript: MANETranscript,
        uta: UTADatabase,
        gene_normalizer: GeneQueryHandler,
        vrs: VRSRepresentation,
        hgvs_dup_del_mode: HGVSDupDelMode
    ) -> None:
        """Initialize the DelIns validator.

        :param seqrepo_access: Access to SeqRepo data
        :param transcript_mappings: Access to transcript mappings
        :param gene_symbol: Gene symbol tokenizer
        :param mane_transcript: Access MANE Transcript information
        :param uta: Access to UTA queries
        :param gene_normalizer: Access to gene-normalizer
        :param vrs: Class for creating VRS objects
        """
        self.transcript_mappings = transcript_mappings
        self.seqrepo_access = seqrepo_access
        self._gene_matcher = gene_symbol
        self.uta = uta
        self.genomic_base = GenomicBase(self.seqrepo_access, self.uta)
        self.mane_transcript = mane_transcript
        self.gene_normalizer = gene_normalizer
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

    async def get_grch38_data(self, classification, errors, ac):
        """

        :param ac: Accession
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

        return {
            "ac": new_ac,
            "pos0": pos0,
            "pos1": pos1
        }

    async def get_grch38_data_ambiguous(
        self, classification, errors
    ) -> Dict:
        old_ac = classification.ac
        pos0, pos1, pos2, pos3, new_ac = None, None, None, None, None
        if classification.ambiguous_type == AmbiguousType.AMBIGUOUS_1:
            grch38_pos0_pos1 = await self.mane_transcript.g_to_grch38(
                old_ac, classification.pos0, classification.pos1
            )
            pos0, pos1 = grch38_pos0_pos1["pos"]
            ac_pos0_pos1 = grch38_pos0_pos1["ac"]

            grch38_pos2_pos3 = await self.mane_transcript.g_to_grch38(
                old_ac, classification.pos2, classification.pos3
            )
            pos2, pos3 = grch38_pos2_pos3["pos"]
            ac_pos2_pos3 = grch38_pos2_pos3["ac"]

            if ac_pos0_pos1 != ac_pos2_pos3:
                errors.append(
                    f"{ac_pos0_pos1} does not equal {ac_pos2_pos3} when lifting over "
                    f"to GRCh38"
                )
            else:
                new_ac = ac_pos0_pos1
        elif classification.ambiguous_type in {AmbiguousType.AMBIGUOUS_2,
                                               AmbiguousType.AMBIGUOUS_5}:
            grch38 = await self.mane_transcript.g_to_grch38(
                old_ac, classification.pos1, classification.pos2
            )
            pos1, pos2 = grch38["pos"]
            new_ac = grch38["ac"]
        elif classification.ambiguous_type == AmbiguousType.AMBIGUOUS_7:
            grch38 = await self.mane_transcript.g_to_grch38(
                old_ac, classification.pos0, classification.pos2
            )
            pos0, pos2, = grch38["pos"]
            new_ac = grch38["ac"]

        if not new_ac:
            errors.append(f"Unable to find a GRCh38 accession for: {old_ac}")

        return {
            "ac": new_ac,
            "pos0": pos0,
            "pos1": pos1,
            "pos2": pos2,
            "pos3": pos3
        }

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
