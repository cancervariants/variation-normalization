from typing import Optional

from biocommons.seqrepo import SeqRepo

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, generic_dna

import hgvs.dataproviders.uta

class SeqRepoAccess:
    #/Users/acoffman/git/varlex/varlexapp/data/seqrepo/latest
    def __init__(self, seqrepo_data_path):
        self.seq_repo_client = SeqRepo(seqrepo_data_path)
        self.hdp = hgvs.dataproviders.uta.connect()


    def transcripts_for_gene_symbol(self, symbol):
        tx = self.hdp.get_tx_for_gene(symbol)

        #vs splign (refseq)? also can filter out coding vs non coding?
        return [a for a in tx if a['alt_aln_method'] == 'genebuild']


    #TODO:definitely will want to cache this
    def protein_at_position(self, transcript, pos: int) -> Optional[str]:
        bases = self.seq_repo_client.fetch(
                transcript['tx_ac'],
                transcript['cds_start_i'],
                transcript['cds_end_i'])

        coding_dna = Seq(bases, generic_dna)
        proteins = coding_dna.translate()
        if len(proteins) >= pos -1:
            return proteins[pos - 1]
        else:
            return None
