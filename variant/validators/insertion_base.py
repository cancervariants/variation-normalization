"""The module for Insertion Validation."""
from variant.schemas.token_response_schema import Token
from variant.validators.validator import Validator
import logging

logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)


class InsertionBase(Validator):
    """The Insertion Validator Base class."""

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

                if allele:
                    self.check_pos_index(t, s, errors)

                self.add_validation_result(
                    allele, valid_alleles, results,
                    classification, s, t, gene_tokens, errors
                )

        # Now add Mane transcripts to results
        self.add_mane_transcript(classification, results, gene_tokens,
                                 mane_transcripts_dict)

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
        """Return a description of the identified variant."""
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
