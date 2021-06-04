"""The module for Genomic Insertion Validation."""
from variation.validators.insertion_base import InsertionBase
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import GenomicInsertionToken
from typing import List, Optional
from variation.schemas.token_response_schema import GeneMatchToken
import logging


logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class GenomicInsertion(InsertionBase):
    """The Genomic Insertion Validator class."""

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

    def get_gene_tokens(self, classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_gene_symbol_tokens(classification)

    def variation_name(self):
        """Return the variation name."""
        return 'genomic insertion'

    def is_token_instance(self, t):
        """Check that token is Genomic Insertion."""
        return t.token_type == 'GenomicInsertion'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is
        Genomic Insertion.
        """
        return classification_type == ClassificationType.GENOMIC_INSERTION

    def human_description(self, transcript,
                          token: GenomicInsertionToken) -> str:
        """Return a human description of the identified variation."""
        if token.inserted_sequence2 is not None:
            inserted_sequence = \
                f"{token.inserted_sequence}_{token.inserted_sequence2}"
        else:
            inserted_sequence = f"{token.inserted_sequence}"
        return f"The insertion of nucleotide(s) {inserted_sequence}" \
               f"between nucleotides g.{token.start_pos_flank} and " \
               f"g.{token.end_pos_flank}"
