"""The module for DelIns Validation."""
from abc import abstractmethod
from variant.validators.validator import Validator
import logging

logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)


class DelInsBase(Validator):
    """The DelIns Validator Base class."""

    def get_allele_from_transcript(self, classification, t, s, errors):
        """Return allele from a given transcript.

        :param Classification classification: The classification for input str
        :param Classification s: Classification token
        :param str t: Transcript
        :param list errors: List of errors
        :return: Allele as a dictionary
        """
        allele = None
        if t.startswith('ENST'):
            return allele

        hgvs_expr, _ = self.get_hgvs_expr(classification, t, s, False)
        return self.get_allele_from_hgvs(hgvs_expr, errors)

    @abstractmethod
    def get_hgvs_expr(self, classification, t, s, is_hgvs):
        """Return a HGVS expression.

        :param Classification classification: A classification for a list of
            tokens
        :param str t: Transcript retrieved from transcript mapping
        """
        raise NotImplementedError

    def concise_description(self, transcript, token) -> str:
        """Return a description of the identified variant."""
        if token.start_pos_del is not None and token.end_pos_del is not None:
            position = f"{token.start_pos_del} to {token.end_pos_del}"
        else:
            position = token.start_pos_del

        if token.inserted_sequence1 is not None and \
                token.inserted_sequence2 is not None:
            sequence = f"{token.inserted_sequence1} to " \
                       f"{token.inserted_sequence2}"
        else:
            sequence = token.inserted_sequence1

        return f'{transcript}:{token.reference_sequence}.' \
               f'{position}delins{sequence}'
