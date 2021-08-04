"""The module for Genomic Silent Mutation Validation."""
from typing import Optional, List
from .single_nucleotide_variation_base import SingleNucleotideVariationBase
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import GenomicSilentMutationToken
import logging

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class GenomicSilentMutation(SingleNucleotideVariationBase):
    """The Genomic Silent Mutation Validator class."""

    def get_transcripts(self, gene_tokens, classification, errors)\
            -> Optional[List[str]]:
        """Get transcript accessions for a given classification.

        :param list gene_tokens: A list of gene tokens
        :param Classification classification: A classification for a list of
            tokens
        :param list errors: List of errors
        :return: List of transcript accessions
        """
        return self.get_genomic_transcripts(classification, errors)

    def get_valid_invalid_results(self, classification_tokens, transcripts,
                                  classification, results,
                                  gene_tokens, normalize_endpoint,
                                  mane_data_found, is_identifier) -> None:
        """Add validation result objects to a list of results.

        :param list classification_tokens: A list of classification Tokens
        :param list transcripts: A list of transcript accessions
        :param Classification classification: A classification for a list of
            tokens
        :param list results: Stores validation result objects
        :param list gene_tokens: List of GeneMatchTokens for a classification
        :param bool normalize_endpoint: `True` if normalize endpoint is being
            used. `False` otherwise.
        :param dict mane_data_found: MANE Transcript information found
        :param bool is_identifier: `True` if identifier is given for exact
            location. `False` otherwise.
        """
        self.silent_mutation_valid_invalid_results(
            classification_tokens, transcripts, classification, results,
            gene_tokens, normalize_endpoint, mane_data_found, is_identifier
        )

    def get_gene_tokens(self, classification):
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_gene_symbol_tokens(classification)

    def variation_name(self):
        """Return the variation name."""
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
        """Return a human description of the identified variation."""
        return f'A genomic silent mutation from {token.position} was a ' \
               f'{token.ref_nucleotide} (the nucleotide was not changed)'
