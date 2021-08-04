"""The module for Coding DNA Insertion Validation."""
from variation.validators.insertion_base import InsertionBase
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import CodingDNAInsertionToken
from typing import List, Optional
from variation.schemas.token_response_schema import GeneMatchToken
import logging


logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class CodingDNAInsertion(InsertionBase):
    """The Coding DNA Insertion Validator class."""

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

    def get_gene_tokens(self, classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_coding_dna_gene_symbol_tokens(classification)

    def variation_name(self):
        """Return the variation name."""
        return 'coding dna insertion'

    def is_token_instance(self, t):
        """Check that token is Coding DNA Insertion."""
        return t.token_type == 'CodingDNAInsertion'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is
        Coding DNA Insertion.
        """
        return classification_type == ClassificationType.CODING_DNA_INSERTION

    def human_description(self, transcript,
                          token: CodingDNAInsertionToken) -> str:
        """Return a human description of the identified variation."""
        if token.inserted_sequence2 is not None:
            inserted_sequence = \
                f"{token.inserted_sequence}_{token.inserted_sequence2}"
        else:
            inserted_sequence = f"{token.inserted_sequence}"
        return f"The insertion of nucleotide(s) {inserted_sequence}" \
               f"between nucleotides c.{token.start_pos_flank} and " \
               f"c.{token.end_pos_flank}"
