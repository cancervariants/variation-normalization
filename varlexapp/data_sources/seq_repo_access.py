"""A module for accessing SeqRepo."""
from typing import Optional
from biocommons.seqrepo import SeqRepo


class SeqRepoAccess:
    """The SeqRepoAccess class."""

    def __init__(self, seqrepo_data_path):
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
                return t[pos - 1]
        except KeyError:
            return None

    def translate_alias(self, input_str):
        """Get Ensembl aliases for gene symbols."""
        aliases = list()
        res = self.seq_repo_client.translate_alias(input_str.strip())
        for alias in res:
            if alias['namespace'] == 'Ensembl':
                aliases.append(alias['alias'])
        return aliases
