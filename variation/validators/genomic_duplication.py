"""The module for Genomic Duplication Validation."""
from .validator import Validator
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import \
    TokenType, DuplicationAltType,\
    GenomicDuplicationToken, GenomicDuplicationRangeToken  # noqa: F401
from typing import List, Optional
from variation.schemas.token_response_schema import GeneMatchToken
import logging
from ga4gh.vrs import models


logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class GenomicDuplication(Validator):
    """The Genomic Duplication Validator class."""

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

    def get_valid_invalid_results(self, classification_tokens, transcripts,
                                  classification, results, gene_tokens,
                                  normalize_endpoint, mane_data_found,
                                  is_identifier) -> None:
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
        for s in classification_tokens:
            for t in transcripts:
                errors = list()
                t = self.get_accession(t, classification)

                start = None
                end = None
                variation = None

                if s.token_type == TokenType.GENOMIC_DUPLICATION:
                    start = s.start_pos1_dup
                    if s.start_pos2_dup is None:
                        # Format: #dup
                        end = s.start_pos1_dup
                    else:
                        # Format: #_#dup
                        end = s.start_pos2_dup

                    variation = self.to_vrs_allele(
                        t, start, end, s.reference_sequence,
                        s.alt_type, errors)
                elif s.token_type == TokenType.GENOMIC_DUPLICATION_RANGE:
                    # TODO: Check if ranges should be CNVs or Alleles
                    if s.alt_type != DuplicationAltType.UNCERTAIN_DUPLICATION:
                        ival = models.SequenceInterval(
                            start=models.DefiniteRange(
                                min=s.start_pos1_dup - 1,
                                max=s.start_pos2_dup - 1),
                            end=models.DefiniteRange(
                                min=s.end_pos1_dup + 1,
                                max=s.end_pos2_dup + 1)
                        )
                        allele = self.to_vrs_allele_ranges(
                            s, t, s.reference_sequence, s.alt_type,
                            errors, ival)
                    else:
                        if s.start_pos1_dup == '?' and s.end_pos2_dup == '?':
                            start = s.start_pos2_dup
                            end = s.end_pos1_dup

                            ival = models.SequenceInterval(
                                start=models.IndefiniteRange(
                                    value=start - 1,
                                    comparator="<="
                                ),
                                end=models.IndefiniteRange(
                                    value=end,
                                    comparator=">="
                                )
                            )
                            allele = self.to_vrs_allele_ranges(
                                s, t, s.reference_sequence, s.alt_type,
                                errors, ival)
                        elif s.start_pos1_dup == '?' and \
                                s.start_pos2_dup != '?' and \
                                s.end_pos1_dup != '?' and \
                                s.end_pos2_dup is None:
                            # format: (?_#)_#
                            start = s.start_pos2_dup
                            end = s.end_pos1_dup

                            ival = models.SequenceInterval(
                                start=models.IndefiniteRange(
                                    value=start - 1,
                                    comparator="<="
                                ),
                                end=models.Number(value=end)
                            )
                            allele = self.to_vrs_allele_ranges(
                                s, t, s.reference_sequence, s.alt_type,
                                errors, ival
                            )
                        elif s.start_pos1_dup != '?' and \
                                s.start_pos2_dup is None and \
                                s.end_pos1_dup != '?' and\
                                s.end_pos2_dup == '?':
                            # format: #_(#_?)
                            start = s.end_pos1_dup
                            end = s.end_pos1_dup

                            ival = models.SequenceInterval(
                                start=models.Number(value=start),
                                end=models.IndefiniteRange(
                                    value=end,
                                    comparator=">="
                                ),

                            )
                            allele = self.to_vrs_allele_ranges(
                                s, t, s.reference_sequence, s.alt_type,
                                errors, ival
                            )
                        else:
                            allele = None
                    if allele is None:
                        errors.append("Unable to get allele")
                        return None
                    variation = self.to_vrs_cnv(t, allele, 'dup')
                    if not variation:
                        errors.append(f"Unable to get CNV for {t}")
                else:
                    errors.append(f"Token type not supported: {s.token_type}")

                if not errors:
                    if not gene_tokens:
                        if s.token_type == TokenType.GENOMIC_DUPLICATION_RANGE\
                                and s.alt_type != DuplicationAltType.UNCERTAIN_DUPLICATION:  # noqa: E501
                            # TODO
                            mane = None
                        else:
                            # If no gene tokens, get GRCh38
                            grch38 = self.mane_transcript.g_to_grch38(
                                t, start, end)

                            if grch38:
                                mane = dict(
                                    gene=None,
                                    refseq=grch38['ac'] if grch38['ac'].startswith('NC') else None,  # noqa: E501
                                    ensembl=grch38['ac'] if grch38['ac'].startswith('ENSG') else None,  # noqa: E501
                                    pos=grch38['pos'],
                                    strand=None,
                                    status='GRCh38'
                                )
                            else:
                                mane = None

                        self.add_mane_data(
                            mane, mane_data_found, s.reference_sequence,
                            s.alt_type, s, gene_tokens
                        )
                self.add_validation_result(
                    variation, valid_alleles, results,
                    classification, s, t, gene_tokens, errors
                )

                if is_identifier:
                    break

        self.add_mane_to_validation_results(
            mane_data_found, valid_alleles, results,
            classification, gene_tokens
        )

    def get_gene_tokens(self, classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_gene_symbol_tokens(classification)

    def variation_name(self):
        """Return the variation name."""
        return 'genomic duplication'

    def is_token_instance(self, t):
        """Check that token is an instance of Genomic Duplication."""
        return t.token_type in ['GenomicDuplication',
                                'GenomicDuplicationRange']

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is
        Genomic Duplication.
        """
        return classification_type == \
            ClassificationType.GENOMIC_DUPLICATION

    def human_description(self, transcript,
                          token) -> str:
        """Return a human description of the identified variation."""
        if token.token_type == 'GenomicDuplication':
            descr = "A Genomic Deletion "
        else:
            # Genomic Duplication Range
            descr = "A Genomic Deletion Range "
        return descr

    def concise_description(self, transcript, token):
        """Return a concise description of the identified variation."""
        return "TODO"
