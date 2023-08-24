"""The module for Genomic Reference Agree Validation."""
from typing import List, Optional, Union

from variation.schemas.classification_response_schema import (
    Classification,
    ClassificationType,
    GenomicReferenceAgreeClassification,
    Nomenclature,
)
from variation.schemas.service_schema import ClinVarAssembly
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import Validator


class GenomicReferenceAgree(Validator):
    """The Genomic Reference Agree Validator class."""

    async def get_valid_invalid_results(
        self, classification: GenomicReferenceAgreeClassification, accessions: List[str]
    ) -> List[ValidationResult]:
        """Get list of validation results for a given classification and accessions

        :param classification: A classification for a list of tokens
        :param accessions: A list of accessions for a classification
        :return: List of validation results containing invalid and valid results
        """
        validation_results = []

        for alt_ac in accessions:
            errors = []

            if classification.nomenclature == Nomenclature.GNOMAD_VCF:
                token = classification.matching_tokens[0]
                ref = token.ref
                start_pos = token.pos
                end_pos = token.pos + len(ref)

                invalid_ref_msg = self.validate_reference_sequence(
                    alt_ac, start_pos, end_pos, ref
                )
                if invalid_ref_msg:
                    errors.append(invalid_ref_msg)
            else:
                invalid_ac_pos_msg = self.validate_ac_and_pos(
                    alt_ac, classification.pos
                )
                if invalid_ac_pos_msg:
                    errors.append(invalid_ac_pos_msg)

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
        """Return whether or not the classification type is genomic reference agree"""
        return classification_type == ClassificationType.GENOMIC_REFERENCE_AGREE

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
