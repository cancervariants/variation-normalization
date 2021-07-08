"""The module for Insertion Validation."""
from variation.schemas.token_response_schema import Token
from variation.validators.validator import Validator
import logging

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class InsertionBase(Validator):
    """The Insertion Validator Base class."""

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

        # TODO: Fix MANE for genomic

        for s in classification_tokens:
            for t in transcripts:
                errors = list()

                allele, t, hgvs_expr, is_ensembl = \
                    self.get_allele_with_context(classification, t, s, errors)

                if allele:
                    self.check_pos_index(t, s, errors)
                else:
                    errors.append("Unable to find allele")

                if not errors:
                    mane = self.mane_transcript.get_mane_transcript(
                        t, s.start_pos_flank, s.end_pos_flank,
                        s.reference_sequence,
                        normalize_endpoint=normalize_endpoint
                    )
                    if mane:
                        if s.reference_sequence == 'g':
                            if not gene_tokens and mane['gene']:
                                gene_tokens.append(
                                    self._gene_matcher.match(mane['gene'])
                                )

                            if mane['status'] != 'GRCh38':
                                s.molecule_context = 'transcript'
                                s.reference_sequence = 'c'

                        prefix = f"{mane['refseq']}:" \
                                 f"{s.reference_sequence.lower()}."
                        position = f"{mane['pos'][0]}_{mane['pos'][1]}"
                        if s.inserted_sequence2 is not None:
                            inserted_seq = f"{s.inserted_sequence}_" \
                                           f"{s.inserted_sequence2}"
                        else:
                            inserted_seq = f"{s.inserted_sequence}"
                        mane_hgvs_expr =\
                            f"{prefix}{position}ins{inserted_seq}"
                        self.add_mane_data(mane_hgvs_expr, mane, mane_data,
                                           s)

                self.add_validation_result(
                    allele, valid_alleles, results,
                    classification, s, t, gene_tokens, errors
                )

                if is_hgvs:
                    break

        self.add_mane_to_validation_results(
            mane_data, valid_alleles, results, classification, gene_tokens
        )

    def get_hgvs_expr(self, classification, t, s, is_hgvs):
        """Return a HGVS expression.

        :param Classification classification: A classification for a list of
            tokens
        :param str t: Transcript retrieved from transcript mapping
        :param Token s: The classification token
        :param bool is_hgvs: Whether or not classification is HGVS token
        :return: hgvs expression
        """
        if not is_hgvs:
            prefix = f"{t}:{s.reference_sequence.lower()}."
            position = f"{s.start_pos_flank}_{s.end_pos_flank}"
            if s.inserted_sequence2 is not None:
                inserted_sequence = \
                    f"{s.inserted_sequence}_{s.inserted_sequence2}"
            else:
                inserted_sequence = f"{s.inserted_sequence}"

            hgvs_expr = f"{prefix}{position}ins{inserted_sequence}"
        else:
            hgvs_token = [t for t in classification.all_tokens if
                          isinstance(t, Token) and t.token_type
                          in ['HGVS', 'ReferenceSequence']][0]
            hgvs_expr = hgvs_token.input_string
        return hgvs_expr

    def concise_description(self, transcript, token) -> str:
        """Return a HGVS description of the identified variation.

        :param str transcript: Transcript accession
        :param Token token: Classification token
        :return: HGVS expression
        """
        position = f"{token.start_pos_flank}_{token.end_pos_flank}"
        if token.inserted_sequence2 is not None:
            inserted_sequence = \
                f"{token.inserted_sequence}_{token.inserted_sequence2}"
        else:
            inserted_sequence = f"{token.inserted_sequence}"

        return f'{transcript}:{token.reference_sequence}.' \
               f'{position}ins{inserted_sequence}'

    def check_pos_index(self, t, s, errors):
        """Check that position exists on transcript.

        :param str t: Transcript accession
        :param Token s: Classification token
        :param list errors: List of errors
        """
        len_of_seq = self.seqrepo_access.len_of_sequence(t)
        is_len_lte_start = len_of_seq <= int(s.start_pos_flank)
        is_len_lte_end = len_of_seq <= int(s.end_pos_flank)

        if is_len_lte_end or is_len_lte_start:
            errors.append('Sequence index error')
