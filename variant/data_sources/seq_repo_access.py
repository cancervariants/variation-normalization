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

    def protein_at_position(self, transcript: str, pos: int) -> Optional[str]:
        """Get the protein at a position."""
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

    def aliases(self, input_str):
        """Get aliases for gene symbols."""
        try:
            return self.seq_repo_client.translate_alias(input_str.strip())
        except KeyError:
            return []
