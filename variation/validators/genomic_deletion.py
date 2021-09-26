"""The module for Genomic Deletion Validation."""
from variation.validators.deletion_base import DeletionBase
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import GenomicDeletionToken
from typing import List, Optional
from variation.schemas.token_response_schema import GeneMatchToken
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

    def get_valid_invalid_results(self, classification_tokens, transcripts,
                                  classification, results, gene_tokens,
                                  normalize_endpoint, mane_data_found,
                                  is_identifier, hgvs_dup_del_mode)\
            -> None:
        """Add validation result objects to a list of results.

        :param list classification_tokens: A list of classification Tokens
        :param list transcripts: A list of transcript accessions
        :param Classification classification: A classification for a list of
            tokens
        :param list results: Stores validation result objects
        :param list gene_tokens: List of GeneMatchTokens for a classification
        :param bool normalize_endpoint: `True` if normalize endpoint is being
            used. `False` otherwise.
        :param dict mane_data_found: MANE Transcript information found
        :param bool is_identifier: `True` if identifier is given for exact
            location. `False` otherwise.
        :param str hgvs_dup_del_mode: Must be: `default`, `cnv`,
            `repeated_seq_expr`, `literal_seq_expr`.
            This parameter determines how to represent HGVS dup/del expressions
            as VRS objects.
        """
        valid_alleles = list()
        for s in classification_tokens:
            for t in transcripts:
                errors = list()
                t = self.get_accession(t, classification)
                start, end = s.start_pos_del, None
                if s.end_pos_del is None:
                    end = start
                else:
                    end = s.end_pos_del

                allele = self.to_vrs_allele(t, s.start_pos_del, s.end_pos_del,
                                            s.reference_sequence, s.alt_type,
                                            errors)

                variation = self.hgvs_dup_del_mode.interpret_variation(
                    t, s.alt_type, allele, errors, hgvs_dup_del_mode,
                    pos=(start, end)
                )

                if not errors:
                    self.check_reference_sequence(t, s, errors)

                if not errors:
                    self._get_normalize_variation(
                        gene_tokens, s, t, errors, hgvs_dup_del_mode,
                        mane_data_found, start, end, normalize_endpoint)

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

    def _get_normalize_variation(self, gene_tokens, s, t, errors,
                                 hgvs_dup_del_mode, mane_data_found, start,
                                 end, normalize_endpoint) -> None:
        """Get variation that will be returned in normalize endpoint.

        :param list gene_tokens: List of gene tokens
        :param Token s: Classification token
        :param str t: Accession
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Mode to use for
            interpreting HGVS duplications and deletions
        :param dict mane_data_found: MANE Transcript data found for given query
        :param int start: Start pos change
        :param int end: End pos change
        :param bool normalize_endpoint: `True` if normalize endpoint is being
            used. `False` otherwise.
        """
        if not gene_tokens:
            grch38 = self.mane_transcript.g_to_grch38(t, start, end)

            if grch38:
                for pos in [grch38['pos'][0], grch38['pos'][1]]:
                    self._check_index(grch38['ac'], pos, errors)

                if not errors:
                    allele = self.to_vrs_allele(
                        grch38['ac'], grch38['pos'][0],
                        grch38['pos'][1], s.reference_sequence,
                        s.alt_type, errors
                    )
                    grch38_variation = \
                        self.hgvs_dup_del_mode.interpret_variation(
                            grch38['ac'], s.alt_type, allele,
                            errors, hgvs_dup_del_mode
                        )

                    if grch38_variation:
                        self._add_dict_to_mane_data(
                            grch38['ac'], s, grch38_variation,
                            mane_data_found, 'GRCh38'
                        )
        else:
            # TODO
            pass
            # mane = self.mane_transcript.get_mane_transcript(
            #     t, s.start_pos_del, s.end_pos_del,
            #     s.reference_sequence,
            #     gene=gene_tokens[0].token if gene_tokens else None,
            #     normalize_endpoint=normalize_endpoint
            # )
            #
            # self.add_mane_data(
            #     mane, mane_data_found, s.reference_sequence,
            #     s.alt_type, s, gene_tokens
            # )

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
