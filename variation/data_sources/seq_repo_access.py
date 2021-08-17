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
        self.cache = dict()

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
            return self.seq_repo_client.fetch(transcript, start=start - 1,
                                              end=end)
        except TypeError:
            try:
                start = int(start)
                end = int(end)
                return self.seq_repo_client.fetch(transcript,
                                                  start=start - 1, end=end)
            except ValueError as e:
                logger.warning(e)
                return None
        except KeyError:
            logger.warning(f"Accession {transcript} not found in SeqRepo")
            return None
        except ValueError as e:
            logger.warning(f"{transcript}: {e}")
            return None

    def aliases(self, input_str) -> List[str]:
        """Get aliases for a given input."""
        try:
            return self.seq_repo_client.translate_alias(input_str.strip())
        except KeyError:
            logger.warning(f"SeqRepo could not translate alias: {input_str}")
            return []
