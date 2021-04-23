"""A module for accessing SeqRepo."""
from typing import Optional
from biocommons.seqrepo import SeqRepo
from variant import SEQREPO_DATA_PATH


class SeqRepoAccess:
    """The SeqRepoAccess class."""

    def __init__(self, seqrepo_data_path=SEQREPO_DATA_PATH):
        """Initialize the SeqRepoAccess class.

        :param str seqrepo_data_path: The path to the seqrepo directory.
        """
        self.seq_repo_client = SeqRepo(seqrepo_data_path)

    def sequence_at_position(self, transcript: str, pos: int) -> Optional[str]:
        """Return sequence at a given transcript position.

        :param str transcript: Transcript
        :param int pos: The position to search on
        :return: A sequence (protein or nucleotide)
        """
        # why does this not exist sometimes?
        try:
            t = self.seq_repo_client.fetch(transcript)
            if len(t) < pos - 1:
                return None
            else:
                try:
                    return t[pos - 1]
                except IndexError:
                    return None
        except KeyError:
            return None

    def len_of_sequence(self, transcript: str):
        """Return the length of a transcript's sequence.

        :param str transcript: Transcript to find sequence length of
        """
        return len(self.seq_repo_client.fetch(transcript))

    def aliases(self, input_str):
        """Get aliases for a given input."""
        try:
            return self.seq_repo_client.translate_alias(input_str.strip())
        except KeyError:
            return []
