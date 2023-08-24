"""The module for Genomic Deletion Ambiguous Validation."""
from typing import List, Optional, Union

from variation.schemas.classification_response_schema import (
    AmbiguousType,
    Classification,
    ClassificationType,
    GenomicDeletionAmbiguousClassification,
    Nomenclature,
)
from variation.schemas.service_schema import ClinVarAssembly
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import Validator


class GenomicDeletionAmbiguous(Validator):
    """The Genomic Deletion Ambiguous Validator class."""

    async def get_valid_invalid_results(
        self,
        classification: GenomicDeletionAmbiguousClassification,
        accessions: List[str],
    ) -> List[ValidationResult]:
        """Get list of validation results for a given classification and accessions

        :param classification: A classification for a list of tokens
        :param accessions: A list of accessions for a classification
        :return: List of validation results containing invalid and valid results
        """
        # Validate ambiguous type and positions
        invalid_classification_msg = self.validate_ambiguous_classification(
            classification
        )
        if invalid_classification_msg:
            return [
                ValidationResult(
                    accession=None,
                    classification=classification,
                    is_valid=False,
                    errors=[invalid_classification_msg],
                )
            ]

        validation_results = []

        for alt_ac in accessions:
            errors = []

            if classification.ambiguous_type == AmbiguousType.AMBIGUOUS_1:
                start_pos = classification.pos0
                end_pos = classification.pos3
            elif classification.ambiguous_type == AmbiguousType.AMBIGUOUS_2:
                start_pos = classification.pos1
                end_pos = classification.pos2
            elif classification.ambiguous_type == AmbiguousType.AMBIGUOUS_5:
                start_pos = classification.pos1
                end_pos = classification.pos2
            elif classification.ambiguous_type == AmbiguousType.AMBIGUOUS_7:
                start_pos = classification.pos0
                end_pos = classification.pos2
            else:
                start_pos = None
                end_pos = None
                errors.append(
                    f"ambiguous type not supported: {classification.ambiguous_type}"
                )

            if start_pos is not None and end_pos is not None:
                invalid_ac_pos = self.validate_ac_and_pos(
                    alt_ac, start_pos, end_pos=end_pos
                )
                if invalid_ac_pos:
                    errors.append(invalid_ac_pos)

            if not errors and classification.gene_token:
                invalid_gene_pos_msg = await self._validate_gene_pos(
                    classification.gene_token.matched_value,
                    alt_ac,
                    classification.pos0,
                    classification.pos1,
                    pos2=classification.pos2,
                    pos3=classification.pos3,
                )
                if invalid_gene_pos_msg:
                    errors.append(invalid_gene_pos_msg)

            validation_results.append(
                ValidationResult(
                    accession=alt_ac,
                    classification=classification,
                    is_valid=not errors,
                    errors=errors,
                )
            )

        return validation_results

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is genomic deletion
        ambiguous
        """
        return classification_type == ClassificationType.GENOMIC_DELETION_AMBIGUOUS

    async def get_accessions(
        self,
        classification: Classification,
        errors: List,
        input_assembly: Optional[
            Union[ClinVarAssembly.GRCH37, ClinVarAssembly.GRCH38]
        ] = None,
    ) -> List[str]:
        """Get accessions for a given classification.
        If `classification.nomenclature == Nomenclature.HGVS`, will return the accession
        in the HGVS expression.
        Else, will get all accessions associated to the gene

        :param classification: The classification for list of tokens
        :param errors: List of errors
        :param input_assembly: Assembly used for initial input query. Only used when
            initial query is using genomic free text or gnomad vcf format
        :return: List of accessions
        """
        if classification.nomenclature == Nomenclature.HGVS:
            accessions = [classification.ac]
        else:
            accessions = await self.get_genomic_accessions(
                classification, errors, input_assembly=input_assembly
            )
        return accessions
