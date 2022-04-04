"""The module for Coding DNA Insertion Validation."""
from typing import List, Optional
import logging

from variation.validators.insertion_base import InsertionBase
from variation.schemas.classification_response_schema import \
    Classification, ClassificationType
from variation.schemas.token_response_schema import Token
from variation.schemas.token_response_schema import GeneMatchToken

logger = logging.getLogger("variation")
logger.setLevel(logging.DEBUG)


class CodingDNAInsertion(InsertionBase):
    """The Coding DNA Insertion Validator class."""

    async def get_transcripts(self, gene_tokens: List, classification: Classification,
                              errors: List) -> Optional[List[str]]:
        """Get transcript accessions for a given classification.

        :param List gene_tokens: A list of gene tokens
        :param Classification classification: A classification for a list of
            tokens
        :param List errors: List of errors
        :return: List of transcript accessions
        """
        return self.get_coding_dna_transcripts(gene_tokens, errors)

    def get_gene_tokens(self, classification: Classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_coding_dna_gene_symbol_tokens(classification)

    def variation_name(self) -> str:
        """Return the variation name."""
        return "coding dna insertion"

    def is_token_instance(self, t: Token) -> bool:
        """Check that token is Coding DNA Insertion."""
        return t.token_type == "CodingDNAInsertion"

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is
        Coding DNA Insertion.
        """
        return classification_type == ClassificationType.CODING_DNA_INSERTION
