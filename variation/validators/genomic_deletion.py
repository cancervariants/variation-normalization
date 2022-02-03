"""The module for Genomic Deletion Validation."""
from variation.validators.duplication_deletion_base import\
    DuplicationDeletionBase
from variation.schemas.classification_response_schema import \
    ClassificationType, Classification
from variation.schemas.token_response_schema import GenomicDeletionToken, \
    Token, SequenceOntology
from typing import List, Optional, Dict
from variation.schemas.token_response_schema import GeneMatchToken
import logging
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum


logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class GenomicDeletion(DuplicationDeletionBase):
    """The Genomic Deletion Validator class."""

    def get_transcripts(self, gene_tokens: List,
                        classification: Classification,
                        errors: List) -> Optional[List[str]]:
        """Get transcript accessions for a given classification.

        :param list gene_tokens: A list of gene tokens
        :param Classification classification: A classification for a list of
            tokens
        :param list errors: List of errors
        :return: List of transcript accessions
        """
        return self.get_genomic_transcripts(classification, errors)

    def get_valid_invalid_results(
            self, classification_tokens: List, transcripts: List,
            classification: Classification, results: List, gene_tokens: List,
            normalize_endpoint: bool, mane_data_found: Dict,
            is_identifier: bool, hgvs_dup_del_mode: HGVSDupDelModeEnum
    ) -> None:
        """Add validation result objects to a list of results.

        :param List classification_tokens: A list of classification Tokens
        :param List transcripts: A list of transcript accessions
        :param Classification classification: A classification for a list of
            tokens
        :param List results: Stores validation result objects
        :param List gene_tokens: List of GeneMatchTokens for a classification
        :param bool normalize_endpoint: `True` if normalize endpoint is being
            used. `False` otherwise.
        :param Dict mane_data_found: MANE Transcript information found
        :param bool is_identifier: `True` if identifier is given for exact
            location. `False` otherwise.
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Must be: `default`, `cnv`,
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

                # Validate pos
                if gene_tokens:
                    self.validate_gene_or_accession_pos(
                        t, [start, end], errors, gene=gene_tokens[0].token)
                else:
                    self.check_reference_sequence(t, s, errors)

                if not errors:
                    allele = self.vrs.to_vrs_allele(
                        t, s.start_pos_del, s.end_pos_del,
                        s.reference_sequence, s.alt_type, errors)

                    variation = self.hgvs_dup_del_mode.interpret_variation(
                        t, s.alt_type, allele, errors, hgvs_dup_del_mode,
                        pos=(start, end)
                    )
                else:
                    variation = None

                if not errors and normalize_endpoint:
                    self._get_normalize_variation(
                        gene_tokens, s, t, errors, hgvs_dup_del_mode,
                        mane_data_found, start, end)

                self.add_validation_result(
                    variation, valid_alleles, results,
                    classification, s, t, gene_tokens, errors
                )

                if is_identifier:
                    break

        if normalize_endpoint:
            self.add_mane_to_validation_results(
                mane_data_found, valid_alleles, results,
                classification, gene_tokens
            )

    def _get_normalize_variation(
            self, gene_tokens: List, s: Token, t: str, errors: List,
            hgvs_dup_del_mode: HGVSDupDelModeEnum, mane_data_found: Dict,
            start: int, end: int) -> None:
        """Get variation that will be returned in normalize endpoint.

        :param List gene_tokens: List of gene tokens
        :param Token s: Classification token
        :param str t: Accession
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Mode to use for
            interpreting HGVS duplications and deletions
        :param Dict mane_data_found: MANE Transcript data found for given query
        :param int start: Start pos change
        :param int end: End pos change
        """
        if gene_tokens:
            self.add_normalized_genomic_dup_del(
                s, t, s.start_pos_del, s.end_pos_del, gene_tokens[0].token,
                SequenceOntology.DELETION, errors, hgvs_dup_del_mode,
                mane_data_found)
        else:
            # No gene provided, then use GRCh38 assesmbly
            if not self._is_grch38_assembly(t):
                grch38 = self.mane_transcript.g_to_grch38(t, start, end)
            else:
                grch38 = dict(ac=t, pos=(start, end))

            if grch38:
                self.validate_gene_or_accession_pos(
                    grch38['ac'], [grch38['pos'][0], grch38['pos'][1]], errors)
                self.add_grch38_to_mane_data(
                    t, s, errors, grch38, mane_data_found, hgvs_dup_del_mode,
                    use_vrs_allele_range=False)

    def get_gene_tokens(
            self, classification: Classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_gene_symbol_tokens(classification)

    def variation_name(self) -> str:
        """Return the variation name."""
        return 'genomic deletion'

    def is_token_instance(self, t: Token):
        """Check that token is Genomic Deletion.

        :param Token t: Classification token
        """
        return t.token_type == 'GenomicDeletion'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is
        Genomic DelIns.

        :param ClassificationType classification_type: Classification type
        :return: `True` if classification type matches, `False` otherwise
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
