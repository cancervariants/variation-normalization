"""The module for Deletion Validation."""
from abc import abstractmethod
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

    def concise_description(self, transcript, token) -> str:
        """Return a description of the identified variant."""
        position = f"{token.start_pos_del}"
        if token.end_pos_del is not None:
            position += f"_{token.end_pos_del}"

        descr = f"{transcript}:{token.reference_sequence}.{position}del"
        if token.deleted_sequence:
            descr += f"{token.deleted_sequence}"
        return descr
