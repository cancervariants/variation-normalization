"""The module for Genomic Silent Mutation Validation."""
from typing import Optional, List

from .single_nucleotide_variant_base import SingleNucleotideVariantBase
from variant.schemas.classification_response_schema import \
    ClassificationType
from variant.schemas.token_response_schema import GenomicSilentMutationToken
import logging
from variant.schemas.token_response_schema import Token

logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)

# TODO: Find gene from NC accession (in event of no mane transcripts)


class GenomicSilentMutation(SingleNucleotideVariantBase):
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
