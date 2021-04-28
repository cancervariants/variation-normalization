"""The module for Genomic Silent Mutation Validation."""
from .single_nucleotide_variant_base import SingleNucleotideVariantBase
from variant.schemas.classification_response_schema import \
    ClassificationType
from variant.schemas.token_response_schema import GenomicSilentMutationToken
from typing import List
from variant.schemas.classification_response_schema import Classification
from variant.schemas.validation_response_schema import ValidationResult
import logging
from variant.schemas.token_response_schema import Token
from .genomic_base import GenomicBase

logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)

# TODO: Find gene from NC accession (in event of no mane transcripts)


class GenomicSilentMutation(SingleNucleotideVariantBase):
    """The Genomic Silent Mutation Validator class."""

    def validate(self, classification: Classification) \
            -> List[ValidationResult]:
        """Validate a given classification.

        :param Classification classification: A classification for a list of
            tokens
        :return: A list of validation results
        """
        results = list()
        errors = list()

        classification_tokens = self.get_classification_tokens(classification)
        gene_tokens = self.get_gene_tokens(classification)

        if gene_tokens and len(gene_tokens) > 1:
            errors.append('More than one gene symbol found for a single'
                          f' {self.variant_name()}')

        if len(classification.non_matching_tokens) > 0:
            errors.append(f"Non matching tokens found for "
                          f"{self.variant_name()}.")

        genomic_base = GenomicBase(self.dp)
        nc_accessions = genomic_base.get_nc_accessions(classification)
        if not nc_accessions:
            errors.append('Could not find NC_ accession for '
                          f'{self.variant_name()}')

        if len(errors) > 0:
            return [self.get_validation_result(
                classification, False, 0, None,
                '', '', errors, gene_tokens)]

        self.get_valid_invalid_results(classification_tokens, nc_accessions,
                                       classification, results, gene_tokens)
        return results

    def get_hgvs_expr(self, classification, t, s, is_hgvs):
        """Get HGVS expression."""
        if is_hgvs:
            hgvs_token = [t for t in classification.all_tokens if
                          isinstance(t, Token) and t.token_type == 'HGVS'][0]
            t = hgvs_token.input_string.split(':')[0]

        hgvs_expr = f"{t}:{s.reference_sequence.lower()}.{s.position}="

        return hgvs_expr, None

    def get_valid_invalid_results(self, classification_tokens, transcripts,
                                  classification, results,
                                  gene_tokens) -> None:
        """Add validation result objects to a list of results.

        :param list classification_tokens: A list of Tokens
        :param list transcripts: A list of transcript strings
        :param Classification classification: A classification for a list of
            tokens
        :param list results: A list to store validation result objects
        :param list gene_tokens: List of GeneMatchTokens
        """
        self.silent_mutation_valid_invalid_results(
            classification_tokens, transcripts, classification, results,
            gene_tokens
        )

    def get_gene_tokens(self, classification):
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_gene_symbol_tokens(classification)

    def variant_name(self):
        """Return the variant name."""
        return 'genomic silent mutation'

    def is_token_instance(self, t):
        """Check that token is Genomic Silent Mutation."""
        return t.token_type == 'GenomicSilentMutation'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is amino acid
        substitution.
        """
        return classification_type == ClassificationType.GENOMIC_SILENT_MUTATION  # noqa: E501

    def human_description(self, transcript,
                          token: GenomicSilentMutationToken) -> str:
        """Return a human description of the identified variant."""
        return f'A genomic silent mutation from {token.position} was a ' \
               f'{token.ref_nucleotide} (the nucleotide was not changed)'
