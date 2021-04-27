"""The module for DelIns Validation."""
from abc import abstractmethod
from variant.validators.validator import Validator
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify
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

    def get_vrs_allele(self, sequence_id, s) -> dict:
        """Return VRS Allele object.

        :param str sequence_id: Sequence containing the sequence to be located
        :param Token s: A Classification token
        :return: A VRS Allele object as a dictionary
        """
        seq_location = models.SequenceLocation(
            sequence_id=sequence_id,
            interval=models.SimpleInterval(
                start=int(s.start_pos_del) - 1,
                end=int(s.end_pos_del)
            )
        )

        state = models.SequenceState(sequence=s.reference_sequence.upper())
        allele = models.Allele(location=seq_location, state=state)
        allele['_id'] = ga4gh_identify(allele)
        return allele.as_dict()

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
