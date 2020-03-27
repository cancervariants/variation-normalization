from typing import Optional

from biocommons.seqrepo import SeqRepo

class SeqRepoAccess:
    #/Users/acoffman/git/varlex/varlexapp/data/seqrepo/latest
    def __init__(self, seqrepo_data_path):
        self.seq_repo_client = SeqRepo(seqrepo_data_path)


    def protein_at_position(self, transcript: str, pos: int) -> Optional[str]:
        #why does this not exist sometimes?
        try:
            t = self.seq_repo_client.fetch(transcript)
            if len(t) < pos -1:
                return None
            else:
                return t[pos - 1]
        except KeyError:
            return None

