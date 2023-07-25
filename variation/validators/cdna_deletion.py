"""The module for cDNA Deletion Validation."""
from typing import List

from variation.schemas.classification_response_schema import (
    ClassificationType, Classification, CdnaDeletionClassification, Nomenclature
)
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.validator import Validator


class CdnaDeletion(Validator):
    """The cDNA Deletion Validator class."""

    async def get_valid_invalid_results(
        self, classification: CdnaDeletionClassification, accessions: List[str]
    ) -> List[ValidationResult]:
        if classification.pos1 and classification.pos0 >= classification.pos1:
            return [ValidationResult(
                accession=None,
                classification=classification,
                is_valid=False,
                errors=[(
                    "Positions deleted should contain two different positions and "
                    "should be listed from 5' to 3'")]
            )]

        validation_results = []

        for c_ac in accessions:
            errors = []
            cds_start, cds_start_err_msg = await self.get_cds_start(c_ac)

            if cds_start_err_msg:
                errors.append(cds_start_err_msg)
            else:
                if classification.nomenclature in {Nomenclature.FREE_TEXT,
                                                   Nomenclature.HGVS}:
                    # # validate deleted sequence
                    # HGVS deleted sequence includes start and end
                    if classification.deleted_sequence:
                        invalid_del_seq_msg = self.validate_reference_sequence(
                            c_ac, cds_start + classification.pos0,
                            cds_start + classification.pos1 + 1,
                            classification.deleted_sequence
                        )

                        if invalid_del_seq_msg:
                            errors.append(invalid_del_seq_msg)
                    else:
                        # Validate accession and positions
                        invalid_ac_pos_msg = self.validate_ac_and_pos(
                            c_ac, cds_start + classification.pos0,
                            end_pos=cds_start + classification.pos1 if classification.pos1 else None  # noqa: E501
                        )
                        if invalid_ac_pos_msg:
                            errors.append(invalid_ac_pos_msg)

            validation_results.append(
                ValidationResult(
                    accession=c_ac,
                    classification=classification,
                    cds_start=cds_start,
                    is_valid=not errors,
                    errors=errors
                )
            )

        return validation_results

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is cdna deletion."""
        return classification_type == ClassificationType.CDNA_DELETION

    async def get_accessions(
        self, classification: Classification, errors: List
    ) -> List[str]:
        """Get accessions for a given classification.
        If `classification.nomenclature == Nomenclature.HGVS`, will return the accession
        in the HGVS expression.
        Else, will get all accessions associated to the gene

        :param classification: The classification for list of tokens
        :param errors: List of errors
        :return: List of accessions
        """
        if classification.nomenclature == Nomenclature.HGVS:
            accessions = [classification.ac]
        else:
            accessions = self.get_cdna_accessions(
                classification.gene_token, errors
            )
        return accessions
