"""Generate genomic accessions cache.
This is a temporary solution for speeding up genomic accession queries.
"""
from variation import APP_ROOT
import json


class GenomicAccessionCache:
    """Class for generating genomic accession cache."""

    def __init__(self, uta, seqrepo) -> None:
        """Initialize GenomicAccessionCache instance.

        :param UTA uta: UTA instance
        :param SeqRepoAccess seqrepo: SeqRepoAccess instance
        """
        self.uta = uta
        self.seqrepo = seqrepo
        self._cache_path = APP_ROOT / 'data' / 'genomic_accessions.json'
        if not self._cache_path.exists():
            self.create_genomic_accessions_json()
        self.cache = self._load_cache()

    def create_genomic_accessions_json(self) -> None:
        """Create JSON for genomic accessions."""
        query = (
            """
            SELECT DISTINCT alt_ac
            FROM uta_20210129.tx_exon_aln_v
            WHERE alt_ac LIKE 'NC_%'
            AND alt_aln_method = 'splign'
            ORDER BY alt_ac
            """
        )
        self.uta.cursor.execute(query)
        genomic_acs = self.uta.cursor.fetchall()
        genomic_acs = [item for sublist in genomic_acs for item in sublist]

        genomic_acs_cache = dict()
        for genomic_ac in genomic_acs:
            genomic_acs_cache[genomic_ac] = self.seqrepo.seq_repo_client.fetch(
                genomic_ac)

        with open(self._cache_path, "w+") as f:
            json.dump(genomic_acs_cache, f)

    def _load_cache(self) -> dict:
        """Load cache containing genomic accessions sequences.

        :return: Genomic accession cache as a dict
        """
        with open(self._cache_path, 'r') as f:
            return json.load(f)
