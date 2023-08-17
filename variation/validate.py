"""Module for Validation."""
from typing import List, Optional, Union

from cool_seq_tool.data_sources import SeqRepoAccess, TranscriptMappings, UTADatabase
from gene.query import QueryHandler as GeneQueryHandler

from variation.schemas.classification_response_schema import Classification
from variation.schemas.service_schema import ClinVarAssembly
from variation.schemas.validation_response_schema import ValidationSummary
from variation.validators import (
    Amplification,
    CdnaDeletion,
    CdnaDelIns,
    CdnaInsertion,
    CdnaReferenceAgree,
    CdnaSubstitution,
    GenomicDeletion,
    GenomicDeletionAmbiguous,
    GenomicDelIns,
    GenomicDuplication,
    GenomicDuplicationAmbiguous,
    GenomicInsertion,
    GenomicReferenceAgree,
    GenomicSubstitution,
    ProteinDeletion,
    ProteinDelIns,
    ProteinInsertion,
    ProteinReferenceAgree,
    ProteinStopGain,
    ProteinSubstitution,
)
from variation.validators.validator import Validator


class Validate:
    """The validation class."""

    def __init__(
        self,
        seqrepo_access: SeqRepoAccess,
        transcript_mappings: TranscriptMappings,
        uta: UTADatabase,
        gene_normalizer: GeneQueryHandler,
    ) -> None:
        """Initialize the validate class. Will create an instance variable,
        `validators`, which is a list of Validators for supported variation types.

        :param seqrepo_access: Access to SeqRepo data
        :param transcript_mappings: Access to transcript mappings
        :param uta: Access to UTA queries
        :param gene_normalizer: Access to gene-normalizer
        """
        params = [seqrepo_access, transcript_mappings, uta, gene_normalizer]
        self.validators: List[Validator] = [
            ProteinSubstitution(*params),
            CdnaSubstitution(*params),
            GenomicSubstitution(*params),
            ProteinStopGain(*params),
            ProteinReferenceAgree(*params),
            CdnaReferenceAgree(*params),
            GenomicReferenceAgree(*params),
            ProteinDelIns(*params),
            CdnaDelIns(*params),
            GenomicDelIns(*params),
            ProteinDeletion(*params),
            CdnaDeletion(*params),
            GenomicDeletion(*params),
            GenomicDeletionAmbiguous(*params),
            ProteinInsertion(*params),
            CdnaInsertion(*params),
            GenomicInsertion(*params),
            GenomicDuplication(*params),
            GenomicDuplicationAmbiguous(*params),
            Amplification(*params),
        ]

    async def perform(
        self,
        classification: Classification,
        input_assembly: Optional[
            Union[ClinVarAssembly.GRCH37, ClinVarAssembly.GRCH38]
        ] = None,
    ) -> ValidationSummary:
        """Get validation summary containing invalid and valid results for a
        classification

        :param classification: A classification for a list of tokens
        :param input_assembly: Assembly used for `q`. Only used when `q` is using
            genomic free text or gnomad vcf format
        :return: Validation summary for classification containing valid and invalid
            results
        """
        valid_possibilities = []
        invalid_possibilities = []

        found_valid_result = False
        invalid_classification = None

        for validator in self.validators:
            if validator.validates_classification_type(
                classification.classification_type
            ):
                validation_results = await validator.validate(
                    classification, input_assembly=input_assembly
                )
                for validation_result in validation_results:
                    if validation_result.is_valid:
                        found_valid_result = True
                        valid_possibilities.append(validation_result)
                    else:
                        invalid_possibilities.append(validation_result)
                        invalid_classification = (
                            classification.classification_type.value
                        )

            if found_valid_result:
                break

        if not found_valid_result:
            warnings = [
                f"Unable to find valid result for classification: {invalid_classification}"  # noqa: E501
            ]
        else:
            warnings = []

        return ValidationSummary(
            valid_results=valid_possibilities,
            invalid_results=invalid_possibilities,
            warnings=warnings,
        )
