"""The module for Coding DNA Substitution Validation."""
from .single_nucleotide_variation_base import SingleNucleotideVariationBase
from variation.schemas.classification_response_schema import \
    ClassificationType, Classification
from variation.schemas.token_response_schema import\
    CodingDNASilentMutationToken
from typing import List, Optional, Dict
from variation.schemas.token_response_schema import GeneMatchToken
import logging
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum


logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class CodingDNASilentMutation(SingleNucleotideVariationBase):
    """The Coding DNA Silent Mutation Validator class."""

    def get_transcripts(self, gene_tokens, classification, errors)\
            -> Optional[List[str]]:
        """Get transcript accessions for a given classification.

        :param list gene_tokens: A list of gene tokens
        :param Classification classification: A classification for a list of
            tokens
        :param list errors: List of errors
        :return: List of transcript accessions
        """
        return self.get_coding_dna_transcripts(gene_tokens, errors)

    def get_valid_invalid_results(
            self, classification_tokens: List, transcripts: List,
            classification: Classification, results: List, gene_tokens: List,
            normalize_endpoint: bool, mane_data_found: Dict,
            is_identifier: bool, hgvs_dup_del_mode: HGVSDupDelModeEnum
    ) -> None:
        """Add validation result objects to a list of results.

        :param List classification_tokens: A list of classification Tokens
        :param List transcripts: A list of transcript accessions
        :param Classification classification: A classification for a list of
            tokens
        :param List results: Stores validation result objects
        :param List gene_tokens: List of GeneMatchTokens for a classification
        :param bool normalize_endpoint: `True` if normalize endpoint is being
            used. `False` otherwise.
        :param Dict mane_data_found: MANE Transcript information found
        :param bool is_identifier: `True` if identifier is given for exact
            location. `False` otherwise.
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Must be: `default`, `cnv`,
            `repeated_seq_expr`, `literal_seq_expr`.
            This parameter determines how to represent HGVS dup/del expressions
            as VRS objects.
        """
        self.silent_mutation_valid_invalid_results(
            classification_tokens, transcripts, classification, results,
            gene_tokens, normalize_endpoint, mane_data_found, is_identifier
        )

    def get_gene_tokens(self, classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_coding_dna_gene_symbol_tokens(classification)

    def variation_name(self):
        """Return the variation name."""
        return 'coding dna silent mutation'

    def is_token_instance(self, t):
        """Check that token is Coding DNA Silent Mutation."""
        return t.token_type == 'CodingDNASilentMutation'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is coding dna silent
        mutation.
        """
        return classification_type == ClassificationType.CODING_DNA_SILENT_MUTATION  # noqa: E501

    def human_description(self, transcript,
                          token: CodingDNASilentMutationToken) -> str:
        """Return a human description of the identified variation."""
        return f'A coding DNA silent mutation from {token.position} ' \
               f'was a {token.ref_nucleotide} (the nucleotide was not changed)'
