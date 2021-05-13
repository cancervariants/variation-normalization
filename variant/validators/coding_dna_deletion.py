"""The module for Coding DNA Deletion Validation."""
from variant.validators.deletion_base import DeletionBase
from variant.schemas.classification_response_schema import \
    ClassificationType
from variant.schemas.token_response_schema import CodingDNADeletionToken
from variant.schemas.validation_response_schema import LookupType
from typing import List
from variant.schemas.classification_response_schema import Classification
from variant.schemas.token_response_schema import GeneMatchToken
from variant.schemas.validation_response_schema import ValidationResult
import logging


logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)


class CodingDNADeletion(DeletionBase):
    """The Coding DNA Deletion Validator class."""

    def validate(self, classification: Classification) \
            -> List[ValidationResult]:
        """Validate a given classification.

        :param Classification classification: A classification for a list of
            tokens
        :return: A list of validation results
        """
        results = list()
        errors = list()

        classification_tokens = self.get_classification_tokens(classification)
        gene_tokens = self.get_gene_tokens(classification)

        if len(classification.non_matching_tokens) > 0:
            errors.append(f"Non matching tokens found for "
                          f"{self.variant_name()}.")

        if len(gene_tokens) == 0:
            errors.append(f'No gene tokens for a {self.variant_name()}.')

        if len(gene_tokens) > 1:
            errors.append('More than one gene symbol found for a single'
                          f' {self.variant_name()}')

        if len(errors) > 0:
            return [self.get_validation_result(
                classification, False, 0, None,
                '', '', errors, gene_tokens)]

        transcripts = self.transcript_mappings.coding_dna_transcripts(
            gene_tokens[0].token, LookupType.GENE_SYMBOL)

        if not transcripts:
            errors.append(f'No transcripts found for gene symbol '
                          f'{gene_tokens[0].token}')
            return [self.get_validation_result(
                classification, False, 0, None,
                '', '', errors, gene_tokens)]

        self.get_valid_invalid_results(classification_tokens, transcripts,
                                       classification, results, gene_tokens)
        return results

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

                if 'HGVS' in classification.matching_tokens:
                    hgvs_expr, is_ensembl_transcript = \
                        self.get_hgvs_expr(classification, t, s, True)
                    allele = self.get_allele_from_hgvs(hgvs_expr, errors)
                    if allele:
                        t = hgvs_expr.split(':')[0]
                    else:
                        errors = list()
                        hgvs_expr, is_ensembl_transcript = \
                            self.get_hgvs_expr(classification, t, s, False)
                        allele = self.get_allele_from_hgvs(hgvs_expr, errors)
                else:
                    hgvs_expr, is_ensembl_transcript = self.get_hgvs_expr(
                        classification, t, s, False
                    )
                    allele = self.get_allele_from_hgvs(hgvs_expr, errors)

                mane_transcripts_dict[hgvs_expr] = {
                    'classification_token': s,
                    'transcript_token': t,
                    'is_ensembl_transcript': is_ensembl_transcript
                }

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

        # Now add Mane transcripts to results
        self.add_mane_transcript(classification, results, gene_tokens,
                                 mane_transcripts_dict)

    def get_gene_tokens(self, classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_coding_dna_gene_symbol_tokens(classification)

    def variant_name(self):
        """Return the variant name."""
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
        """Return a human description of the identified variant."""
        if token.start_pos_del is not None and token.end_pos_del is not None:
            position = f"{token.start_pos_del} to {token.end_pos_del}"
        else:
            position = token.start_pos_del

        descr = "A Coding DNA "
        if token.deleted_sequence:
            descr += f"{token.deleted_sequence} "
        descr += f"Deletion from {position} on transcript {transcript}"
        return descr
