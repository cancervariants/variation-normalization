"""The module for Genomic Duplication Ambiguous Validation."""
from typing import List, Optional

from variation.schemas.classification_response_schema import (
    AmbiguousType, Classification, ClassificationType, GenomicDuplicationAmbiguousClassification,
    Nomenclature
)
from variation.schemas.token_response_schema import GeneToken
from variation.schemas.validation_response_schema import ValidationResult
from .validator import Validator


class GenomicDuplicationAmbiguous(Validator):
    """The Genomic Duplication Ambiguous Validator class."""

    async def get_valid_invalid_results(
        self, classification: GenomicDuplicationAmbiguousClassification,
        transcripts: List[str]
    ) -> List[ValidationResult]:
        if classification.ambiguous_type == AmbiguousType.AMBIGUOUS_1:
            if classification.pos3 <= classification.pos2 <= classification.pos1 <= classification.pos0:  # noqa: E501
                return [ValidationResult(
                    accession=None,
                    classification=classification,
                    is_valid=False,
                    errors=[(
                        "Positions duplicated should contain two different positions and "
                        "should be listed from 5' to 3'")]
                )]
        elif classification.ambiguous_type in {AmbiguousType.AMBIGUOUS_2,
                                               AmbiguousType.AMBIGUOUS_5}:
            if classification.pos2 <= classification.pos1:
                return [ValidationResult(
                    accession=None,
                    classification=classification,
                    is_valid=False,
                    errors=[(
                        "Positions duplicated should contain two different positions and "
                        "should be listed from 5' to 3'")]
                )]
        elif classification.ambiguous_type == AmbiguousType.AMBIGUOUS_7:
            if classification.pos2 <= classification.pos0:
                return [ValidationResult(
                    accession=None,
                    classification=classification,
                    is_valid=False,
                    errors=[(
                        "Positions duplicated should contain two different positions and "
                        "should be listed from 5' to 3'")]
                )]
        else:
            return [ValidationResult(
                accession=None,
                classification=classification,
                is_valid=False,
                errors=[(f"{classification.ambiguous_type} is not yet supported")]
            )]

        validation_results = []
        # _validate_gene_pos?

        for ac in transcripts:
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
                    ac, start_pos, end_pos=end_pos
                )
                if invalid_ac_pos:
                    errors.append(invalid_ac_pos)

            validation_results.append(
                ValidationResult(
                    accession=ac,
                    classification=classification,
                    is_valid=not errors,
                    errors=errors
                )
            )

        return validation_results

    def variation_name(self) -> str:
        """Return the variation name."""
        return "genomic duplication ambiguous"

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is genomic duplication
        ambiguous
        """
        return classification_type == ClassificationType.GENOMIC_DUPLICATION_AMBIGUOUS

    async def get_transcripts(
        self, classification: Classification, errors: List
    ) -> List[str]:
        """Get transcript accessions for a given classification.

        :param Classification classification: A classification for a list of
            tokens
        :param List errors: List of errors
        :return: List of transcript accessions
        """
        if classification.nomenclature == Nomenclature.HGVS:
            transcripts = [classification.ac]
        else:
            transcripts = await self.get_genomic_transcripts(
                classification, errors
            )
        return transcripts
