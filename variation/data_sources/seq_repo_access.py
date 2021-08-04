"""A module for accessing SeqRepo."""
from typing import Optional, List
from biocommons.seqrepo import SeqRepo
from variation import SEQREPO_DATA_PATH
import logging


logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class SeqRepoAccess:
    """The SeqRepoAccess class."""

    def __init__(self, seqrepo_data_path=SEQREPO_DATA_PATH):
        """Initialize the SeqRepoAccess class.

        :param str seqrepo_data_path: The path to the seqrepo directory.
        """
        self.seq_repo_client = SeqRepo(seqrepo_data_path)

    def sequence_at_position(self, transcript: str, pos: int) -> Optional[str]:
        """Return sequence at a position for a given transcript.

        :param str transcript: Transcript accession
        :param int pos: The position to search on
        :return: A sequence (protein or nucleotide)
        """
        try:
            t = self.seq_repo_client.fetch(transcript)
            if len(t) < pos - 1:
                logger.warning(f"Position {pos} exceeds {transcript} length")
                return None
            else:
                try:
                    return t[pos - 1]
                except IndexError:
                    logger.warning(f"Position {pos} not found in transcript "
                                   f"{transcript}")
                    return None
        except KeyError:
            logger.warning(f"Accession {transcript} not found in SeqRepo")
            return None

    def get_sequence(self, transcript: str, start: int, end: int)\
            -> Optional[str]:
        """Return sequence from start and end position for a transcript.

        :param str transcript: Transcript accession
        :param int start: Start position
        :param int end: End position
        """
        try:
            sequence = self.seq_repo_client.fetch(transcript)
            len_of_sequence = len(sequence)
            if len_of_sequence < start - 1:
                logger.warning(f"Position {start-1} exceeds"
                               f" {transcript} length")
                return None
            elif len_of_sequence < end:
                logger.warning(f"Position {end} exceeds"
                               f" {transcript} length")
                return None
            else:
                try:
                    return sequence[start - 1:end]
                except IndexError:
                    logger.warning(f"Positions {start -1} to {end} not found "
                                   f"in transcript {transcript}")
                    return None
        except KeyError:
            logger.warning(f"Accession {transcript} not found in SeqRepo")
            return None

    def len_of_sequence(self, transcript: str) -> int:
        """Return the length of a transcript's sequence.

        :param str transcript: Transcript to find sequence length of
        :return: Length of transcript
        """
        try:
            return len(self.seq_repo_client.fetch(transcript))
        except KeyError:
            logger.warning(f"Accession {transcript} not found in SeqRepo")
            return 0

    def aliases(self, input_str) -> List[str]:
        """Get aliases for a given input."""
        try:
            return self.seq_repo_client.translate_alias(input_str.strip())
        except KeyError:
            logger.warning(f"SeqRepo could not translate alias: {input_str}")
            return []
