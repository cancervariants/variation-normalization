"""Module for translation."""
from abc import abstractmethod, ABC
from typing import Dict, Optional

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
from variation.schemas.classification_response_schema import ClassificationType


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
        vrs: VRSRepresentation
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
    ) -> Optional[Dict]:
        """TODO:"""
