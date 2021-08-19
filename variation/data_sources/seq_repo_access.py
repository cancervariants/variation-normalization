"""A module for accessing SeqRepo."""
from typing import Optional, List
from biocommons.seqrepo import SeqRepo
from variation import SEQREPO_DATA_PATH
from os import environ
import logging


logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class SeqRepoAccess:
    """The SeqRepoAccess class."""

    def __init__(self, seqrepo_data_path=SEQREPO_DATA_PATH):
        """Initialize the SeqRepoAccess class.

        :param str seqrepo_data_path: The path to the seqrepo directory.
        """
        environ['SEQREPO_LRU_CACHE_MAXSIZE'] = "none"
        self.seq_repo_client = SeqRepo(seqrepo_data_path)

    def get_sequence(self, transcript, start, end=None) -> Optional[str]:
        """Return sequence for transcript at given positions.

        :param str transcript: Accession
        :param int start: Start pos change
        :param int end: End pos change
        :return: Sequence
        """
        if end is None:
            end = start

        try:
            sequence = self.seq_repo_client.fetch(transcript, start=start - 1,
                                                  end=end)
            return self.is_valid_index(transcript, end, sequence)
        except TypeError:
            try:
                start = int(start)
                end = int(end)
                sequence = self.seq_repo_client.fetch(transcript,
                                                      start=start - 1, end=end)
                return self.is_valid_index(transcript, end, sequence)
            except ValueError as e:
                logger.warning(e)
                return None
        except KeyError:
            logger.warning(f"Accession {transcript} not found in SeqRepo")
            return None
        except ValueError as e:
            logger.warning(f"{transcript}: {e}")
            return None

    def is_valid_index(self, ac, pos, sequence) -> Optional[str]:
        """Check that index actually exists and return sequence if it does.

        :param str ac: Accession
        :param int pos: End position to check
        :param str sequence: Sequence at pos change
        :return: Sequence at position change
        """
        if self.seq_repo_client.fetch(ac, pos - 1, end=pos):
            return sequence
        else:
            logger.warning(f"Index Error: pos {pos} out of range on {ac}")
            return None

    def aliases(self, input_str) -> List[str]:
        """Get aliases for a given input."""
        try:
            return self.seq_repo_client.translate_alias(input_str.strip())
        except KeyError:
            logger.warning(f"SeqRepo could not translate alias: {input_str}")
            return []
