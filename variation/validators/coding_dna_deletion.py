"""The module for Coding DNA Deletion Validation."""
from variation.validators.deletion_base import DeletionBase
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import CodingDNADeletionToken
from typing import List, Optional
from variation.schemas.token_response_schema import GeneMatchToken
from variation.schemas.token_response_schema import Token
import logging


logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class CodingDNADeletion(DeletionBase):
    """The Coding DNA Deletion Validator class."""

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
        if not is_hgvs:
            prefix = f"{t}:{s.reference_sequence.lower()}.{s.start_pos_del}"
            if s.end_pos_del:
                prefix += f"_{s.end_pos_del}"
            hgvs_expr = f"{prefix}del"
            if s.deleted_sequence:
                hgvs_expr += f"{s.deleted_sequence}"
        else:
            hgvs_token = [t for t in classification.all_tokens if
                          isinstance(t, Token) and t.token_type
                          in ['HGVS', 'ReferenceSequence']][0]
            hgvs_expr = hgvs_token.input_string
        return hgvs_expr

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
        valid_alleles = list()
        mane_transcripts_dict = dict()
        for s in classification_tokens:
            for t in transcripts:
                errors = list()
                allele, t, hgvs_expr, is_ensembl = \
                    self.get_allele_with_context(classification, t, s, errors)

                if hgvs_expr not in mane_transcripts_dict.keys():
                    mane_transcripts_dict[hgvs_expr] = {
                        'classification_token': s,
                        'transcript_token': t,
                        'nucleotide': is_ensembl
                    }

                if not allele:
                    errors.append("Unable to find allele.")
                else:
                    ref_sequence = self.get_reference_sequence(t, s, errors)

                    if ref_sequence and s.deleted_sequence:
                        self.check_reference_sequence(
                            ref_sequence, s.deleted_sequence, errors
                        )

                self.add_validation_result(
                    allele, valid_alleles, results,
                    classification, s, t, gene_tokens, errors
                )

        # Now add Mane transcripts to results
        self.add_mane_transcript(classification, results, gene_tokens,
                                 mane_transcripts_dict)

    def get_gene_tokens(self, classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_coding_dna_gene_symbol_tokens(classification)

    def variation_name(self):
        """Return the variation name."""
        return 'coding dna deletion'

    def is_token_instance(self, t):
        """Check that token is Coding DNA Deletion."""
        return t.token_type == 'CodingDNADeletion'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is
        Coding DNA Deletion.
        """
        return classification_type == ClassificationType.CODING_DNA_DELETION

    def human_description(self, transcript,
                          token: CodingDNADeletionToken) -> str:
        """Return a human description of the identified variation."""
        if token.start_pos_del is not None and token.end_pos_del is not None:
            position = f"{token.start_pos_del} to {token.end_pos_del}"
        else:
            position = token.start_pos_del

        descr = "A Coding DNA "
        if token.deleted_sequence:
            descr += f"{token.deleted_sequence} "
        descr += f"Deletion from {position} on transcript {transcript}"
        return descr
