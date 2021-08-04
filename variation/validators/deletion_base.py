"""The module for Deletion Validation."""
from typing import Optional
from variation.validators.validator import Validator
import logging

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class DeletionBase(Validator):
    """The Deletion Validator Base class."""

    def get_reference_sequence(self, ac, start, end, errors, cds_start=None)\
            -> Optional[str]:
        """Get deleted reference sequence.

        :param str ac: Accession
        :param int start: Start position
        :param int end: End position
        :param list errors: List of errors
        :param int cds_start: Coding start site
        :return: Reference sequence of nucleotides
        """
        if cds_start:
            start += cds_start
            if end is not None:
                end += cds_start

        if start and not end:
            ref_sequence = self.seqrepo_access.sequence_at_position(
                ac, start
            )
        elif start is not None and end is not None:
            ref_sequence = self.seqrepo_access.get_sequence(
                ac, start, end
            )
        else:
            ref_sequence = None

        if not ref_sequence:
            errors.append("Unable to get reference sequence.")
        return ref_sequence

    def check_reference_sequence(self, t, s, errors, cds_start=None) -> bool:
        """Check that reference sequence matches deleted sequence.

        :param str t: Accession
        :param Token s: Classification token
        :param list errors: List of errors
        :param int cds_start: Coding start site
        :return: `True` if ref_sequences matches deleted_sequence.
            `False` otherwise.
        """
        ref_sequence = self.get_reference_sequence(
            t, s.start_pos_del, s.end_pos_del, errors, cds_start=cds_start
        )

        if not errors and ref_sequence and s.deleted_sequence:
            if ref_sequence != s.deleted_sequence:
                errors.append(f"Expected deleted sequence {ref_sequence} "
                              f"but got {s.deleted_sequence}")

    def concise_description(self, transcript, token) -> str:
        """Return a HGVS description of the identified variation.

        :param str transcript: Transcript accession
        :param Token token: Classification token
        :return: HGVS expression
        """
        position = f"{token.start_pos_del}"
        if token.end_pos_del is not None:
            position += f"_{token.end_pos_del}"

        descr = f"{transcript}:{token.reference_sequence}.{position}del"
        if token.deleted_sequence:
            descr += f"{token.deleted_sequence}"
        return descr
