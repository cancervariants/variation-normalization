"""The module for Genomic Insertion Validation."""
from typing import List, Optional
import logging

from variation.validators.insertion_base import InsertionBase
from variation.schemas.classification_response_schema import \
    Classification, ClassificationType
from variation.schemas.token_response_schema import GeneMatchToken, Token


logger = logging.getLogger("variation")
logger.setLevel(logging.DEBUG)


class GenomicInsertion(InsertionBase):
    """The Genomic Insertion Validator class."""

    async def get_transcripts(self, gene_tokens: List, classification: Classification,
                              errors: List) -> Optional[List[str]]:
        """Get transcript accessions for a given classification.

        :param List gene_tokens: A list of gene tokens
        :param Classification classification: A classification for a list of
            tokens
        :param List errors: List of errors
        :return: List of transcript accessions
        """
        transcripts = await self.get_genomic_transcripts(classification, errors)
        return transcripts

    def get_gene_tokens(self, classification: Classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_gene_symbol_tokens(classification)

    def variation_name(self) -> str:
        """Return the variation name."""
        return "genomic insertion"

    def is_token_instance(self, t: Token) -> bool:
        """Check that token is Genomic Insertion."""
        return t.token_type == "GenomicInsertion"

    def validates_classification_type(
            self, classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is
        Genomic Insertion.
        """
        return classification_type == ClassificationType.GENOMIC_INSERTION
