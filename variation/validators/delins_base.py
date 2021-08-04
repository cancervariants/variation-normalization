"""The module for DelIns Validation."""
from variation.validators.validator import Validator
import logging

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class DelInsBase(Validator):
    """The DelIns Validator Base class."""

    def concise_description(self, transcript, token) -> str:
        """Return a HGVS description of the identified variation.

        :param str transcript: Transcript accession
        :param Token token: Classification token
        :return: HGVS expression
        """
        if token.start_pos_del is not None and token.end_pos_del is not None:
            position = f"{token.start_pos_del}_{token.end_pos_del}"
        else:
            position = token.start_pos_del

        if token.inserted_sequence1 is not None and \
                token.inserted_sequence2 is not None:
            sequence = f"{token.inserted_sequence1}_" \
                       f"{token.inserted_sequence2}"
        else:
            sequence = token.inserted_sequence1

        return f'{transcript}:{token.reference_sequence}.' \
               f'{position}delins{sequence}'

    def check_pos_index(self, t, s, errors):
        """Check that position exists on transcript.

        :param str t: Transcript accession
        :param Token s: Classification token
        :param list errors: List of errors
        """
        len_of_seq = self.seqrepo_access.len_of_sequence(t)
        is_len_lte_start = len_of_seq <= int(s.start_pos_del)
        is_len_lte_end = \
            s.end_pos_del and (len_of_seq <= int(s.end_pos_del))

        if is_len_lte_end or is_len_lte_start:
            errors.append('Sequence index error')
