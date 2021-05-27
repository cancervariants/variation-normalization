"""The module for Coding DNA Substitution Validation."""
from .single_nucleotide_variant_base import SingleNucleotideVariantBase
from variant.schemas.classification_response_schema import \
    ClassificationType
from variant.schemas.token_response_schema import CodingDNASilentMutationToken
from typing import List, Optional
from variant.schemas.token_response_schema import GeneMatchToken
from variant.schemas.token_response_schema import Token
import logging


logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)


class CodingDNASilentMutation(SingleNucleotideVariantBase):
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

    def get_hgvs_expr(self, classification, t, s, is_hgvs) -> str:
        """Return HGVS expression

        :param Classification classification: A classification for a list of
            tokens
        :param str t: Transcript retrieved from transcript mapping
        :param Token s: The classification token
         :param bool is_hgvs: Whether or not classification is HGVS token
        :return: hgvs expression
        """
        # Get transcript
        if is_hgvs:
            hgvs_token = [t for t in classification.all_tokens if
                          isinstance(t, Token) and t.token_type == 'HGVS'][0]
            input_string = hgvs_token.input_string
            # TODO: Check if this is ok to do
            if not input_string.startswith('ENST'):
                t = input_string.split(':')[0]
        return f"{t}:{s.reference_sequence.lower()}.{s.position}="

    def get_valid_invalid_results(self, classification_tokens, transcripts,
                                  classification, results, gene_tokens) \
            -> None:
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

    def get_gene_tokens(self, classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_coding_dna_gene_symbol_tokens(classification)

    def variant_name(self):
        """Return the variant name."""
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
        """Return a human description of the identified variant."""
        return f'A coding DNA silent mutation from {token.position} ' \
               f'was a {token.ref_nucleotide} (the nucleotide was not changed)'
