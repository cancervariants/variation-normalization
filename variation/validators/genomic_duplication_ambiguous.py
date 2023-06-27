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
        transcripts: List[str], gene_tokens: List[GeneToken]
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
        # TODO: Validate pos0 and pos1 exist on given accession
        # _validate_gene_pos?

        for ac in transcripts:
            validation_results.append(
                ValidationResult(
                    accession=ac,
                    classification=classification,
                    is_valid=True,
                    errors=[],
                    gene_tokens=gene_tokens
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
        self, gene_tokens: List, classification: Classification, errors: List
    ) -> Optional[List[str]]:
        """Get transcript accessions for a given classification.

        :param List gene_tokens: A list of gene tokens
        :param Classification classification: A classification for a list of
            tokens
        :param List errors: List of errors
        :return: List of transcript accessions
        """
        if classification.nomenclature == Nomenclature.HGVS:
            transcripts = [classification.ac]
        else:
            transcripts = await self.get_genomic_transcripts(
                classification, gene_tokens, errors
            )
        return transcripts

    def get_gene_tokens(self, classification: Classification) -> List:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_gene_symbol_tokens(classification)
