"""The module for Genomic Duplication Validation."""
from .validator import Validator
from variation.schemas.classification_response_schema import \
    ClassificationType
from variation.schemas.token_response_schema import \
    TokenType, DuplicationAltType,\
    GenomicDuplicationToken, GenomicDuplicationRangeToken  # noqa: F401
from typing import List, Optional, Dict, Tuple
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

                result = self._get_variation(s, t, errors, hgvs_dup_del_mode)
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

    def _get_variation(self, s, t, errors, hgvs_dup_del_mode)\
            -> Optional[Dict]:
        """Get variation data.

        :param Token s: Classification token
        :param str t: Accession
        :param list errors: List of errors
        :param HGVSDupDelMode hgvs_dup_del_mode: Mode to use for interpreting
            HGVS duplications and deletions
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
            allele = self.to_vrs_allele(
                t, start, end, s.reference_sequence,
                s.alt_type, errors)
            variation = self._interpret_variation(
                t, s.alt_type, allele, errors, hgvs_dup_del_mode,
                pos=(start, end))
        elif s.token_type == TokenType.GENOMIC_DUPLICATION_RANGE:
            ival, _ = self._get_ival(t, s, errors)

            allele = self.to_vrs_allele_ranges(
                s, t, s.reference_sequence, s.alt_type,
                errors, ival
            )
            if start is not None and end is not None:
                pos = (start, end)
            else:
                pos = None
            variation = self._interpret_variation(
                t, s.alt_type, allele, errors,
                hgvs_dup_del_mode, pos=pos)
        else:
            errors.append(f"Token type not supported: {s.token_type}")

        return {
            'start': start,
            'end': end,
            'variation': variation
        }

    def _add_to_mane_data(self, ac, s, variation, mane_data, status) -> None:
        """Add variation to mane data for normalize endpoint.

        :param str ac: Accession
        :param Token s: Classification token
        :param dict variation: VRS Variation object
        :param dict mane_data: MANE Transcript data found for given query
        :param str status: Status for variation (GRCh38, MANE Select, etc)
        """
        _id = variation['_id']
        key = '_'.join(status.lower().split())

        if _id in mane_data[key].keys():
            mane_data[key][_id]['count'] += 1
        else:
            mane_data[key][_id] = {
                'classification_token': s,
                'accession': ac,
                'count': 1,
                'variation': variation,
                'label': ac  # TODO: Use VRS to translate
            }

    def _get_normalize_variation(self, gene_tokens, s, t, errors,
                                 hgvs_dup_del_mode, mane_data_found, start,
                                 end) -> None:
        """Get variation that will be returned in normalize endpoint.

        :param list gene_tokens: List of gene tokens
        :param Token s: Classification token
        :param str t: Accession
        :param HGVSDupDelMode hgvs_dup_del_mode: Mode to use for interpreting
            HGVS duplications and deletions
        :param dict mane_data_found: MANE Transcript data found for given query
        :param int start: Start pos change
        :param int end: End pos change
        """
        if not gene_tokens:
            if s.token_type == TokenType.GENOMIC_DUPLICATION_RANGE:
                if s.alt_type != DuplicationAltType.UNCERTAIN_DUPLICATION:  # noqa: E501
                    start_grch38 = self.mane_transcript.g_to_grch38(
                        t, s.start_pos1_dup, s.start_pos2_dup
                    )['pos']
                    end_grch38 = self.mane_transcript.g_to_grch38(
                        t, s.end_pos1_dup, s.end_pos2_dup
                    )['pos']

                    for pos in [start_grch38[0], start_grch38[1],
                                end_grch38[0], end_grch38[1]]:
                        self._check_index(t, pos, errors)

                    if not errors:
                        ival = self._get_ival_certain_range(
                            start_grch38[0], start_grch38[1],
                            end_grch38[0], end_grch38[1]
                        )
                        allele = self.to_vrs_allele_ranges(
                            s, t, s.reference_sequence, s.alt_type,
                            errors, ival
                        )
                        mane_variation = self._interpret_variation(
                            t, s.alt_type, allele, errors,
                            hgvs_dup_del_mode)
                        if mane_variation:
                            grch38 = self._grch38_dict(t, None)
                            self.add_mane_data(
                                grch38, mane_data_found, s.reference_sequence,  # noqa: E501
                                s.alt_type, s, gene_tokens, mane_variation=mane_variation  # noqa: E501
                            )
                else:
                    ival, grch38 = self._get_ival(t, s, errors, is_norm=True)

                    if not errors:
                        allele = self.to_vrs_allele_ranges(
                            s, t, s.reference_sequence,
                            s.alt_type, errors, ival)
                        if start is not None and end is not None:
                            pos = (start, end)
                        else:
                            pos = None
                        grch38_variation = self._interpret_variation(
                            t, s.alt_type, allele, errors, hgvs_dup_del_mode
                        )

                        if grch38_variation:
                            self._add_to_mane_data(
                                grch38['ac'], s, grch38_variation,
                                mane_data_found, 'GRCh38'
                            )
            else:
                # If no gene tokens, get GRCh38
                grch38 = self.mane_transcript.g_to_grch38(
                    t, start, end)

                if grch38:
                    for pos in [grch38['pos'][0], grch38['pos'][1]]:  # noqa: E501
                        self._check_index(grch38['ac'], pos, errors)
                    if not errors:
                        allele = self.to_vrs_allele(
                            grch38['ac'], grch38['pos'][0],
                            grch38['pos'][1], s.reference_sequence,
                            s.alt_type, errors
                        )
                        grch38_variation = \
                            self._interpret_variation(
                                grch38['ac'], s.alt_type, allele,
                                errors, hgvs_dup_del_mode
                            )

                        if grch38_variation:
                            self._add_to_mane_data(
                                grch38['ac'], s, grch38_variation,
                                mane_data_found, 'GRCh38'
                            )

    def _grch38_dict(self, ac, pos) -> Dict:
        """Create dict for normalized concepts

        :param str ac: Acession
        :param tuple pos: Position changes
        :return: GRCh38 data
        """
        return dict(
            gene=None,
            refseq=ac if ac.startswith('NC') else None,
            ensembl=ac if ac.startswith('ENSG') else None,
            pos=pos,
            strand=None,
            status='GRCh38'
        )

    def _interpret_variation(self, ac, alt_type, allele, errors,
                             hgvs_dup_del_mode, pos=None) -> Dict:
        """Interpret variation using HGVSDupDelMode

        :param str ac: Accession
        :param str alt_type: Alteration type
        :param dict allele: VRS Allele object
        :param list errors: List of errors
        :param HGVSDupDelMode hgvs_dup_del_mode: Mode to use for interpreting
            HGVS duplications and deletions
        :param tuple pos: Position changes
        :return: VRS Variation object
        """
        variation = None
        if allele is None:
            errors.append("Unable to get Allele")
        else:
            if hgvs_dup_del_mode == 'default':
                variation = self.hgvs_dup_del_mode.default_mode(
                    ac, alt_type, pos, 'dup', allele['location'], allele=allele
                )
            elif hgvs_dup_del_mode == 'cnv':
                variation = self.hgvs_dup_del_mode.cnv_mode(
                    ac, 'dup', allele['location']
                )
            elif hgvs_dup_del_mode == 'repeated_seq_expr':
                variation = self.hgvs_dup_del_mode.repeated_seq_expr_mode(
                    alt_type, allele['location']
                )
            elif hgvs_dup_del_mode == 'literal_seq_expr':
                variation = \
                    self.hgvs_dup_del_mode.literal_seq_expr_mode(allele)
            if not variation:
                errors.append("Unable to get VRS Variation")
        return variation

    def _check_index(self, ac, pos, errors) -> Optional[str]:
        """Check that index actually exists

        :param str ac: Accession
        :param tuple pos: Position changes
        :param list errors: List of errors
        :return: Sequence
        """
        seq = self.seqrepo_access.get_sequence(ac, pos)
        if not seq:
            errors.append(f"Pos {pos} not found on {ac}")
        return seq

    def _get_ival(self, t, s, errors, is_norm=False)\
            -> Optional[Tuple[models.SequenceInterval, Dict]]:
        """Get ival for variations with ranges.

        :param str t: Accession
        :param Token t: Classification token
        :param list errors: List of errors
        :param bool is_norm: `True` if normalize endpoint is being used.
            `False` otherwise.
        :return: Sequence Interval and GRCh38 data if normalize endpoint
            is being used
        """
        ival = None
        start = None
        end = None
        grch38 = None
        if s.alt_type != DuplicationAltType.UNCERTAIN_DUPLICATION:
            # (#_#)_(#_#)
            if is_norm:
                grch38 = self.mane_transcript.g_to_grch38(
                    t, s.start_pos1_dup, s.end_pos2_dup
                )
                if grch38:
                    start, end = grch38['pos']
            else:
                start = s.start_pos1_dup
                end = s.end_pos2_dup

            if start and end:
                ival = self._get_ival_certain_range(
                    s.start_pos1_dup, s.start_pos2_dup,
                    s.end_pos1_dup, s.end_pos2_dup
                )
        else:
            if s.start_pos1_dup == '?' and s.end_pos2_dup == '?':
                # format: (?_#)_(#_?)
                if is_norm:
                    grch38 = self.mane_transcript.g_to_grch38(
                        t, s.start_pos2_dup, s.end_pos1_dup
                    )
                    if grch38:
                        start, end = grch38['pos']
                else:
                    start = s.start_pos2_dup
                    end = s.end_pos1_dup

                if start and end:
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
                    grch38 = self.mane_transcript.g_to_grch38(
                        t, s.start_pos2_dup, s.end_pos1_dup
                    )
                    if grch38:
                        start, end = grch38['pos']
                else:
                    start = s.start_pos2_dup
                    end = s.end_pos1_dup

                if start and end:
                    ival = models.SequenceInterval(
                        start=self._get_start_indef_range(start),  # noqa: E501
                        end=models.Number(value=end)
                    )
            elif s.start_pos1_dup != '?' and \
                    s.start_pos2_dup is None and \
                    s.end_pos1_dup != '?' and \
                    s.end_pos2_dup == '?':
                # format: #_(#_?)
                if is_norm:
                    grch38 = self.mane_transcript.g_to_grch38(
                        t, s.start_pos2_dup, s.end_pos1_dup
                    )
                    if grch38:
                        start, end = grch38['pos']
                else:
                    start = s.end_pos1_dup
                    end = s.end_pos1_dup

                if start and end:
                    ival = models.SequenceInterval(
                        start=models.Number(value=start),
                        end=self._get_end_indef_range(end)
                    )
            else:
                errors.append("Not yet supported")
        return ival, grch38

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
