"""The module for Genomic Uncertain Deletion Validation."""
from variation.validators.duplication_deletion_base import\
    DuplicationDeletionBase
from variation.schemas.classification_response_schema import \
    ClassificationType, Classification
from variation.schemas.token_response_schema import \
    GenomicUncertainDeletionToken, Token
from typing import List, Optional, Dict, Tuple
from variation.schemas.token_response_schema import GeneMatchToken
import logging
from ga4gh.vrs import models
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum


logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class GenomicUncertainDeletion(DuplicationDeletionBase):
    """The Genomic UncertainDeletion Validator class."""

    def get_transcripts(self, gene_tokens: List,
                        classification: Classification,
                        errors: List) -> Optional[List[str]]:
        """Get transcript accessions for a given classification.

        :param List gene_tokens: A list of gene tokens
        :param Classification classification: A classification for a list of
            tokens
        :param List errors: List of errors
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

                result = self._get_variation(s, t, errors, gene_tokens,
                                             hgvs_dup_del_mode)
                variation = result['variation']

                if not errors and normalize_endpoint:
                    self._get_normalize_variation(
                        gene_tokens, s, t, errors, hgvs_dup_del_mode,
                        mane_data_found)

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

    def _get_variation(
            self, s: Token, t: str, errors: List, gene_tokens: List,
            hgvs_dup_del_mode: HGVSDupDelModeEnum) -> Optional[Dict]:
        """Get variation data.

        :param Token s: Classification token
        :param str t: Accession
        :param List errors: List of errors
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Mode to use for
            interpreting HGVS duplications and deletions
        :return: Dictionary containing start/end position changes and variation
        """
        variation, start, end = None, None, None
        ival, grch38 = self._get_ival(t, s, errors, gene_tokens)

        if not errors:
            if grch38:
                t = grch38['ac']

            allele = self.vrs.to_vrs_allele_ranges(
                t, s.reference_sequence, s.alt_type, errors, ival)
            if start is not None and end is not None:
                pos = (start, end)
            else:
                pos = None
            variation = self.hgvs_dup_del_mode.interpret_variation(
                t, s.alt_type, allele, errors,
                hgvs_dup_del_mode, pos=pos)

        return {
            'start': start,
            'end': end,
            'variation': variation
        }

    def _get_normalize_variation(
            self, gene_tokens: List, s: Token, t: str, errors: List,
            hgvs_dup_del_mode: HGVSDupDelModeEnum,
            mane_data_found: Dict) -> None:
        """Get variation that will be returned in normalize endpoint.

        :param List gene_tokens: List of gene tokens
        :param Token s: Classification token
        :param str t: Accession
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Mode to use for
            interpreting HGVS duplications and deletions
        :param dict mane_data_found: MANE Transcript data found for given query
        """
        if not gene_tokens:
            ival, grch38 = self._get_ival(
                t, s, errors, gene_tokens, is_norm=True)
            self.add_grch38_to_mane_data(
                t, s, errors, grch38, mane_data_found, hgvs_dup_del_mode,
                ival=ival)

    def _get_ival(
            self, t: str, s: Token, errors: List, gene_tokens: List,
            is_norm: bool = False
    ) -> Optional[Tuple[models.SequenceInterval, Dict]]:
        """Get ival for variations with ranges.

        :param str t: Accession
        :param Token t: Classification token
        :param List errors: List of errors
        :param bool is_norm: `True` if normalize endpoint is being used.
            `False` otherwise.
        :return: Sequence Interval and GRCh38 data if normalize endpoint
            is being used
        """
        ival = None
        grch38 = None
        gene = gene_tokens[0].token if gene_tokens else None
        if s.start_pos1_del == '?' and s.end_pos2_del == '?':
            # format: (?_#)_(#_?)
            if is_norm:
                t, start, end, _, _, grch38 = self.get_grch38_pos_ac(
                    t, s.start_pos2_del, s.end_pos1_del
                )
            else:
                start = s.start_pos2_del
                end = s.end_pos1_del

            self.validate_gene_or_accession_pos(
                t, [start, end], errors, gene=gene)

            if not errors and start and end:
                ival = models.SequenceInterval(
                    start=self.vrs.get_start_indef_range(start),
                    end=self.vrs.get_end_indef_range(end)
                )
        elif s.start_pos1_del == '?' and \
                s.start_pos2_del != '?' and \
                s.end_pos1_del != '?' and \
                s.end_pos2_del is None:
            # format: (?_#)_#
            if is_norm:
                t, start, end, _, _, grch38 = self.get_grch38_pos_ac(
                    t, s.start_pos2_del, s.end_pos1_del
                )
            else:
                start = s.start_pos2_del
                end = s.end_pos1_del

            self.validate_gene_or_accession_pos(
                t, [start, end], errors, gene=gene
            )

            if not errors and start and end:
                ival = models.SequenceInterval(
                    start=self.vrs.get_start_indef_range(start),  # noqa: E501
                    end=models.Number(value=end)
                )
        elif s.start_pos1_del != '?' and \
                s.start_pos2_del is None and \
                s.end_pos1_del != '?' and \
                s.end_pos2_del == '?':
            # format: #_(#_?)
            if is_norm:
                t, start, end, _, _, grch38 = self.get_grch38_pos_ac(
                    t, s.start_pos1_del, s.end_pos1_del
                )
            else:
                start = s.start_pos1_del
                end = s.end_pos1_del

            start -= 1

            self.validate_gene_or_accession_pos(
                t, [start, end], errors, gene=gene)

            if not errors and start and end:
                ival = models.SequenceInterval(
                    start=models.Number(value=start),
                    end=self.vrs.get_end_indef_range(end)
                )
        else:
            errors.append("Not yet supported")
        return ival, grch38

    def get_gene_tokens(
            self, classification: Classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_gene_symbol_tokens(classification)

    def variation_name(self) -> str:
        """Return the variation name."""
        return 'genomic uncertain deletion'

    def is_token_instance(self, t: Token) -> bool:
        """Check that token is Genomic Uncertain Deletion.

        :param Token t: Classification token
        """
        return t.token_type == 'GenomicUncertainDeletion'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is
        Genomic Uncertain Deletion.

        :param ClassificationType classification_type: Classification type
        :return: `True` if classification type matches, `False` otherwise
        """
        return classification_type == \
            ClassificationType.GENOMIC_UNCERTAIN_DELETION

    def human_description(self, transcript: str,
                          token: GenomicUncertainDeletionToken) -> str:
        """Return a human description of the identified variation.

        :param str transcript: Accession
        :param GenomicUncertainDeletionToken token: Classification token
        """
        descr = f"A Genomic Uncertain Deletion from" \
                f" (?_{token.start_pos2_del}) to {token.end_pos1_del}_? " \
                f"on {transcript}"
        return descr

    def concise_description(self, transcript: str,
                            token: GenomicUncertainDeletionToken) -> str:
        """Return a concise description of the identified variation.

        :param str transcript: Accession
        :param GenomicUncertainDeletionToken token: Classification token
        """
        return f"{transcript}:g.(?_{token.start_pos2_del})_" \
               f"({token.end_pos1_del}_?)"
