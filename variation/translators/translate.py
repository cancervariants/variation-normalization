"""Module for translation."""
from typing import List, Optional, Dict

from gene.query import QueryHandler as GeneQueryHandler
from ga4gh.vrsatile.pydantic.vrs_models import CopyChange
from cool_seq_tool.data_sources import (
    TranscriptMappings, SeqRepoAccess, UTADatabase, MANETranscript
)

from variation.schemas.app_schemas import Endpoint
from variation.schemas.validation_response_schema import ValidationResult
from variation.schemas.normalize_response_schema import (
    HGVSDupDelMode as HGVSDupDelModeEnum
)
from variation.tokenizers import GeneSymbol
from variation.vrs_representation import VRSRepresentation
from .translator import Translator
from .protein_substitution import ProteinSubstitution
from .protein_reference_agree import ProteinReferenceAgree
from .coding_dna_substitution import CdnaSubstitution
from .genomic_substitution import GenomicSubstitution
from .coding_dna_reference_agree import CodingDNAReferenceAgree
from .genomic_reference_agree import GenomicReferenceAgree
from .protein_delins import ProteinDelIns
from .coding_dna_delins import CodingDNADelIns
from .genomic_delins import GenomicDelIns
from .protein_deletion import ProteinDeletion
from .coding_dna_deletion import CdnaDeletion
from .genomic_deletion import GenomicDeletion
from .protein_insertion import ProteinInsertion
from .coding_dna_insertion import CodingDNAInsertion
from .genomic_insertion import GenomicInsertion
from .genomic_uncertain_deletion import GenomicUncertainDeletion
from .genomic_duplication import GenomicDuplication
from .genomic_deletion_range import GenomicDeletionRange
from .amplification import Amplification


class Translate:
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
        """Initialize the translation class."""
        params = [
            seqrepo_access, transcript_mappings, gene_symbol,
            mane_transcript, uta, gene_normalizer, vrs
        ]

        self.all_translators: List[Translator] = [
            ProteinSubstitution(*params),
            CdnaSubstitution(*params),
            GenomicSubstitution(*params),
            ProteinReferenceAgree(*params),
            # CodingDNAReferenceAgree(),
            # GenomicReferenceAgree(),
            ProteinDelIns(*params),
            # CodingDNADelIns(),
            # GenomicDelIns(),
            ProteinDeletion(*params),
            CdnaDeletion(*params),
            # GenomicDeletion(),
            ProteinInsertion(*params),
            # CodingDNAInsertion(),
            GenomicInsertion(*params),
            # GenomicDeletionRange(),
            # GenomicUncertainDeletion(),
            # GenomicDuplication(),
            Amplification(*params)
        ]

    async def perform(
        self,
        validation_result: ValidationResult,  # this is always valid
        warnings: List[str],
        endpoint_name: Optional[Endpoint] = None,
        hgvs_dup_del_mode: HGVSDupDelModeEnum = HGVSDupDelModeEnum.DEFAULT,
        baseline_copies: Optional[int] = None,
        copy_change: Optional[CopyChange] = None,
        do_liftover: bool = False
    ) -> Optional[Dict]:
        """Translate a valid variation query."""
        for translator in self.all_translators:
            if translator.can_translate(
                validation_result.classification.classification_type
            ):
                variation = await translator.translate(
                    validation_result,
                    warnings,
                    endpoint_name=endpoint_name,
                    hgvs_dup_del_mode=hgvs_dup_del_mode,
                    baseline_copies=baseline_copies,
                    copy_change=copy_change,
                    do_liftover=do_liftover
                )
                return variation
        return None
