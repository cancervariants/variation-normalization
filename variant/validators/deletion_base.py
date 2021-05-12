"""The module for Deletion Validation."""
from abc import abstractmethod
from typing import Optional

from variant.validators.validator import Validator
import logging

logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)


class DeletionBase(Validator):
    """The Deletion Validator Base class."""

    @abstractmethod
    def get_hgvs_expr(self, classification, t, s, is_hgvs):
        """Return a HGVS expression.

        :param Classification classification: A classification for a list of
            tokens
        :param str t: Transcript retrieved from transcript mapping
        """
        raise NotImplementedError

    def get_reference_sequence(self, t, s, errors) -> Optional[str]:
        """Get deleted reference sequence.

        :param str t: Transcript
        :param Classification s: Classification for a list of tokens
        :param list errors: List of errors
        :return: Reference sequence of nucleotides
        """
        if s.start_pos_del and not s.end_pos_del:
            ref_sequence = \
                self.seqrepo_access.sequence_at_position(t, s.start_pos_del)
        elif s.start_pos_del and s.end_pos_del:
            ref_sequence = self.seqrepo_access.get_sequence(
                t, s.start_pos_del, s.end_pos_del
            )
        else:
            ref_sequence = None
        if not ref_sequence:
            errors.append("Unable to get reference sequence.")
        return ref_sequence

    @staticmethod
    def check_reference_sequence(ref_sequence,
                                 deleted_sequence, errors) -> bool:
        """Check that reference sequence matches deleted sequence.

        :param str ref_sequence: The reference deleted sequence
        :param str deleted_sequence: The given deleted sequence
        :param list errors: List of errors
        :return: `True` if ref_sequences matches deleted_sequence.
            `False` otherwise.
        """
        if ref_sequence != deleted_sequence:
            errors.append(f"Expected deleted sequence {ref_sequence} "
                          f"but got {deleted_sequence}")

    def concise_description(self, transcript, token) -> str:
        """Return a description of the identified variant."""
        position = f"{token.start_pos_del}"
        if token.end_pos_del is not None:
            position += f"_{token.end_pos_del}"

        descr = f"{transcript}:{token.reference_sequence}.{position}del"
        if token.deleted_sequence:
            descr += f"{token.deleted_sequence}"
        return descr
