"""The module for Genomic Duplication Validation."""
from variation.validators.duplication_deletion_base import\
    DuplicationDeletionBase
from variation.schemas.classification_response_schema import \
    ClassificationType, Classification
from variation.schemas.token_response_schema import \
    TokenType, DuplicationAltType, \
    GenomicDuplicationToken, GenomicDuplicationRangeToken, Token  # noqa: F401
from typing import List, Optional, Dict, Tuple
from variation.schemas.token_response_schema import GeneMatchToken
import logging
from ga4gh.vrs import models
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class GenomicDuplication(DuplicationDeletionBase):
    """The Genomic Duplication Validator class."""

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

                result = self._get_variation(
                    s, t, errors, gene_tokens, hgvs_dup_del_mode,
                    gene=gene_tokens[0].token if gene_tokens else None)
                variation = result['variation']
                start = result['start']
                end = result['end']

                if not errors:
                    self._get_normalize_variation(
                        gene_tokens, s, t, errors, hgvs_dup_del_mode,
                        mane_data_found, start, end)

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

    def _get_variation(self, s: Token, t: str, errors: List, gene_tokens: List,
                       hgvs_dup_del_mode: HGVSDupDelModeEnum,
                       gene: str = None) -> Optional[Dict]:
        """Get variation data.

        :param Token s: Classification token
        :param str t: Accession
        :param List errors: List of errors
        :param HGVSDupDelMode hgvs_dup_del_mode: Mode to use for interpreting
            HGVS duplications and deletions
        :param str gene: Gene symbol token
        :return: Dictionary containing start/end position changes and variation
        """
        variation, start, end = None, None, None
        if s.token_type == TokenType.GENOMIC_DUPLICATION:
            start = s.start_pos1_dup
            if s.start_pos2_dup is None:
                # Format: #dup
                end = s.start_pos1_dup
            else:
                # Format: #_#dup
                end = s.start_pos2_dup

            self.validate_gene_or_accession_pos(
                t, [start, end], errors, gene=gene)

            if not errors:
                allele = self.to_vrs_allele(
                    t, start, end, s.reference_sequence,
                    s.alt_type, errors)
                variation = self.hgvs_dup_del_mode.interpret_variation(
                    t, s.alt_type, allele, errors, hgvs_dup_del_mode,
                    pos=(start, end))
        elif s.token_type == TokenType.GENOMIC_DUPLICATION_RANGE:
            ival, grch38 = self._get_ival(t, s, gene_tokens, errors)

            if not errors:
                if grch38:
                    t = grch38['ac']

                allele = self.to_vrs_allele_ranges(
                    t, s.reference_sequence, s.alt_type, errors, ival)
                if start is not None and end is not None:
                    pos = (start, end)
                else:
                    pos = None
                variation = self.hgvs_dup_del_mode.interpret_variation(
                    t, s.alt_type, allele, errors,
                    hgvs_dup_del_mode, pos=pos)
        else:
            errors.append(f"Token type not supported: {s.token_type}")

        return {
            'start': start,
            'end': end,
            'variation': variation
        }

    def _get_normalize_variation(self, gene_tokens: List, s: Token, t: str,
                                 errors: List,
                                 hgvs_dup_del_mode: HGVSDupDelModeEnum,
                                 mane_data_found: Dict, start: int,
                                 end: int) -> None:
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
        if s.token_type == TokenType.GENOMIC_DUPLICATION_RANGE:
            if s.alt_type != DuplicationAltType.UNCERTAIN_DUPLICATION:
                # (#_#)_(#_#)
                ival, grch38 = self._get_ival(
                    t, s, gene_tokens, errors, is_norm=True)
                self.add_grch38_to_mane_data(
                    t, s, errors, grch38, mane_data_found,
                    hgvs_dup_del_mode, ival=ival
                )
            else:
                ival, grch38 = self._get_ival(t, s, gene_tokens, errors,
                                              is_norm=True)
                self.add_grch38_to_mane_data(
                    t, s, errors, grch38, mane_data_found,
                    hgvs_dup_del_mode, ival=ival)
        else:
            # #dup or #_#dup
            if gene_tokens:
                gene = gene_tokens[0].token

                # Validate position
                self._validate_gene_pos(gene, t, start, end, errors)
                if errors:
                    return

                self.add_normalized_genomic_dup_del(
                    s, t, start, end, gene_tokens[0].token, 'SO:1000035',
                    errors, hgvs_dup_del_mode, mane_data_found)
            else:
                grch38 = self.mane_transcript.g_to_grch38(
                    t, start, end)

                if grch38:
                    self.validate_gene_or_accession_pos(
                        grch38['ac'], [grch38['pos'][0], grch38['pos'][1]],
                        errors)
                    self.add_grch38_to_mane_data(
                        t, s, errors, grch38, mane_data_found,
                        hgvs_dup_del_mode, use_vrs_allele_range=False
                    )

    def _get_ival(
            self, t: str, s: Token, gene_tokens: List, errors: List,
            is_norm: bool = False
    ) -> Optional[Tuple[models.SequenceInterval, Dict]]:
        """Get ival for variations with ranges.

        :param str t: Accession
        :param Token t: Classification token
        :param List gene_tokens: List of gene tokens
        :param List errors: List of errors
        :param bool is_norm: `True` if normalize endpoint is being used.
            `False` otherwise.
        :return: Sequence Interval and GRCh38 data if normalize endpoint
            is being used
        """
        ival, start, end, grch38 = None, None, None, None
        gene = gene_tokens[0].token if gene_tokens else None
        if s.alt_type != DuplicationAltType.UNCERTAIN_DUPLICATION:
            # (#_#)_(#_#)
            if is_norm:
                t, start1, start2, end1, end2, grch38 = self.get_grch38_pos_ac(
                    t, s.start_pos1_dup, s.start_pos2_dup, pos3=s.end_pos1_dup,
                    pos4=s.end_pos2_dup
                )
            else:
                start1 = s.start_pos1_dup
                start2 = s.start_pos2_dup
                end1 = s.end_pos1_dup
                end2 = s.end_pos2_dup

            if start1 and start2 and end1 and end2:
                self.validate_gene_or_accession_pos(
                    t, [start1, start2, end1, end2], errors, gene=gene)

                if not errors:
                    ival = self._get_ival_certain_range(
                        start1, start2, end1, end2)
        else:
            if s.start_pos1_dup == '?' and s.end_pos2_dup == '?':
                # format: (?_#)_(#_?)
                if is_norm:
                    t, start, end, _, _, grch38 = self.get_grch38_pos_ac(
                        t, s.start_pos2_dup, s.end_pos1_dup
                    )
                else:
                    start = s.start_pos2_dup
                    end = s.end_pos1_dup

                # Validate positions
                self.validate_gene_or_accession_pos(
                    t, [start, end], errors, gene=gene)

                if not errors and start and end:
                    ival = models.SequenceInterval(
                        start=self._get_start_indef_range(start),
                        end=self._get_end_indef_range(end)
                    )
            elif s.start_pos1_dup == '?' and \
                    s.start_pos2_dup != '?' and \
                    s.end_pos1_dup != '?' and \
                    s.end_pos2_dup is None:
                # format: (?_#)_#
                if is_norm:
                    t, start, end, _, _, grch38 = self.get_grch38_pos_ac(
                        t, s.start_pos2_dup, s.end_pos1_dup)
                else:
                    start = s.start_pos2_dup
                    end = s.end_pos1_dup

                # Validate positions
                self.validate_gene_or_accession_pos(
                    t, [start, end], errors, gene=gene)

                if not errors and start and end:
                    ival = models.SequenceInterval(
                        start=self._get_start_indef_range(start),
                        end=models.Number(value=end)
                    )
            elif s.start_pos1_dup != '?' and \
                    s.start_pos2_dup is None and \
                    s.end_pos1_dup != '?' and \
                    s.end_pos2_dup == '?':
                # format: #_(#_?)
                if is_norm:
                    t, start, end, _, _, grch38 = self.get_grch38_pos_ac(
                        t, s.start_pos1_dup, s.end_pos1_dup
                    )
                else:
                    start = s.start_pos1_dup
                    end = s.end_pos1_dup
                start -= 1

                self.validate_gene_or_accession_pos(
                    t, [start, end], errors, gene=gene)

                if not errors and start and end:
                    ival = models.SequenceInterval(
                        start=models.Number(value=start),
                        end=self._get_end_indef_range(end)
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

        :param ClassificationType classification_type: Classification type
        :return: `True` if classification type matches, `False` otherwise
        """
        return classification_type == \
            ClassificationType.GENOMIC_DUPLICATION

    def human_description(
            self, transcript: str, token: Token) -> str:
        """Return a human description of the identified variation.

        :param str transcript: Accession
        :param Token token: Classification token
        :return: Human description of variant
        """
        if token.token_type == 'GenomicDuplication':
            descr = "A Genomic Deletion "
        else:
            descr = "A Genomic Deletion Range "
        return descr

    def concise_description(self, transcript: str, token: Token) -> str:
        """Return a concise description of the identified variation.

        :param str transcript: Accession
        :param Token token: Classification token
        """
        return ""
