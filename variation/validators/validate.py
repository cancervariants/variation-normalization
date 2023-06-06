"""Module for Validation."""
from typing import List, Optional

from ga4gh.vrsatile.pydantic.vrs_models import CopyChange
from gene.query import QueryHandler as GeneQueryHandler
from cool_seq_tool.data_sources import TranscriptMappings, SeqRepoAccess, UTADatabase, \
    MANETranscript

from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum
from variation.vrs_representation import VRSRepresentation
from variation.schemas.app_schemas import Endpoint
from variation.schemas.validation_response_schema import ValidationSummary
from variation.schemas.classification_response_schema import Classification
from variation.tokenizers import GeneSymbol
from .protein_substitution import ProteinSubstitution
from .protein_reference_agree import ProteinReferenceAgree
from .coding_dna_substitution import CdnaSubstitution
from .coding_dna_reference_agree import CodingDNAReferenceAgree
from .genomic_reference_agree import GenomicReferenceAgree
from .genomic_substitution import GenomicSubstitution
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


class Validate:
    """The validation class."""

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
        """Initialize the validate class.

        :param seqrepo_access: Access to SeqRepo data
        :param transcript_mappings: Access to transcript mappings
        :param gene_symbol: Gene symbol tokenizer
        :param mane_transcript: Access MANE Transcript information
        :param uta: Access to UTA queries
        :param tlr: Class for translating nomenclatures to and from VRS
        :param gene_normalizer: Access to gene-normalizer
        :param vrs: Class for representing VRS objects
        """
        params = [
            seqrepo_access, transcript_mappings, gene_symbol,
            mane_transcript, uta, gene_normalizer, vrs
        ]
        self.validators = [
            ProteinSubstitution(*params),
            CdnaSubstitution(*params),
            GenomicSubstitution(*params),
            ProteinReferenceAgree(*params),
            # CodingDNAReferenceAgree(*params),
            # GenomicReferenceAgree(*params),
            ProteinDelIns(*params),
            CodingDNADelIns(*params),
            GenomicDelIns(*params),
            ProteinDeletion(*params),
            CdnaDeletion(*params),
            # GenomicDeletion(*params),
            ProteinInsertion(*params),
            # CodingDNAInsertion(*params),
            GenomicInsertion(*params),
            # GenomicDeletionRange(*params),
            # GenomicUncertainDeletion(*params),
            # GenomicDuplication(*params),
            Amplification(*params)
        ]

    async def perform(
        self, classifications: List[Classification], warnings: List = None
    ) -> ValidationSummary:
        valid_possibilities = []
        invalid_possibilities = []
        if not warnings:
            warnings = []

        found_valid_result = False
        invalid_classifications = set()
        for classification in classifications:
            for validator in self.validators:
                if validator.validates_classification_type(
                    classification.classification_type
                ):
                    validation_results = await validator.validate(classification)
                    for validation_result in validation_results:
                        if validation_result.is_valid:
                            found_valid_result = True
                            valid_possibilities.append(validation_result)
                        else:
                            invalid_possibilities.append(validation_result)
                            invalid_classifications.add(
                                classification.classification_type.value
                            )

                if found_valid_result:
                    break

        if not found_valid_result and not warnings:
            warnings.append(f"Unable to find valid result for classifications: "
                            f"{invalid_classifications}")

        return ValidationSummary(
            valid_results=valid_possibilities,
            invalid_results=invalid_possibilities,
            warnings=warnings
        )
