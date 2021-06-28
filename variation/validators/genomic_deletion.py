"""The module for Genomic Deletion Validation."""
from variation.validators.deletion_base import DeletionBase
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import GenomicDeletionToken
from typing import List, Optional
from variation.schemas.token_response_schema import GeneMatchToken
from variation.schemas.token_response_schema import Token
import logging


logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class GenomicDeletion(DeletionBase):
    """The Genomic Deletion Validator class."""

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

    def get_hgvs_expr(self, classification, t, s, is_hgvs) -> str:
        """Return HGVS expression

        :param Classification classification: A classification for a list of
            tokens
        :param str t: Transcript retrieved from transcript mapping
        :param Token s: The classification token
        :param bool is_hgvs: Whether or not classification is HGVS token
        :return: hgvs expression
        """
        if not is_hgvs:
            prefix = f"{t}:{s.reference_sequence.lower()}.{s.start_pos_del}"
            if s.end_pos_del:
                prefix += f"_{s.end_pos_del}"
            hgvs_expr = f"{prefix}del"
            if s.deleted_sequence:
                hgvs_expr += f"{s.deleted_sequence}"
        else:
            hgvs_token = [t for t in classification.all_tokens if
                          isinstance(t, Token) and t.token_type == 'HGVS'][0]
            hgvs_expr = hgvs_token.input_string
        return hgvs_expr

    def get_valid_invalid_results(self, classification_tokens, transcripts,
                                  classification, results, gene_tokens,
                                  normalize_endpoint) -> None:
        """Add validation result objects to a list of results.

        :param list classification_tokens: A list of Tokens
        :param list transcripts: A list of transcript strings
        :param Classification classification: A classification for a list of
            tokens
        :param list results: A list to store validation result objects
        :param list gene_tokens: List of GeneMatchTokens
        :param bool normalize_endpoint: `True` if normalize endpoint is being
            used. `False` otherwise.
        """
        valid_alleles = list()
        if 'HGVS' in classification.matching_tokens:
            is_hgvs = True
        else:
            is_hgvs = False

        mane_data = {
            'mane_select': dict(),
            'mane_plus_clinical': dict(),
            'longest_compatible_remaining': dict()
        }

        for s in classification_tokens:
            for t in transcripts:
                errors = list()
                allele, t, hgvs_expr, is_ensembl = \
                    self.get_allele_with_context(classification, t, s, errors)

                mane = self.mane_transcript.get_mane_transcript(
                    t, s.end_pos_del, s.end_pos_del, s.reference_sequence,
                    normalize_endpoint=normalize_endpoint
                )
                # TODO: Fix MANE when GRCh38 rather than mane
                if mane:
                    if not gene_tokens:
                        gene_tokens.append(
                            self._gene_matcher.match(mane['gene'])
                        )

                    prefix = f"{mane['refseq']}:" \
                             f"c.{mane['pos'][0]}"
                    if s.end_pos_del:
                        prefix += f"_{mane['pos'][1]}"
                    mane_hgvs_expr = f"{prefix}del"
                    if s.deleted_sequence:
                        mane_hgvs_expr += f"{s.deleted_sequence}"
                    self.add_mane_data(mane_hgvs_expr, mane, mane_data, s)

                if allele:
                    ref_sequence = self.get_reference_sequence(t, s, errors)

                    if ref_sequence and s.deleted_sequence:
                        self.check_reference_sequence(
                            ref_sequence, s.deleted_sequence, errors
                        )

                self.add_validation_result(
                    allele, valid_alleles, results,
                    classification, s, t, gene_tokens, errors
                )

                if is_hgvs:
                    break

        self.add_mane_to_validation_results(
            mane_data, valid_alleles, results, classification, gene_tokens
        )

    def get_gene_tokens(self, classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_gene_symbol_tokens(classification)

    def variation_name(self):
        """Return the variation name."""
        return 'genomic deletion'

    def is_token_instance(self, t):
        """Check that token is Genomic Deletion."""
        return t.token_type == 'GenomicDeletion'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is
        Genomic DelIns.
        """
        return classification_type == ClassificationType.GENOMIC_DELETION

    def human_description(self, transcript,
                          token: GenomicDeletionToken) -> str:
        """Return a human description of the identified variation."""
        if token.start_pos_del is not None and token.end_pos_del is not None:
            position = f"{token.start_pos_del} to {token.end_pos_del}"
        else:
            position = token.start_pos_del

        descr = "A Genomic "
        if token.deleted_sequence:
            descr += f"{token.deleted_sequence} "
        descr += f"Deletion from {position} on transcript {transcript}"
        return descr
