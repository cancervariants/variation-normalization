"""Module for Genomic Validation methods."""
from variation.data_sources import UTA
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
import logging


logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class GenomicBase:
    """Genomic Base class for validation methods."""

    def __init__(self, dp: SeqRepoDataProxy, uta: UTA):
        """Initialize the Genomic base class.

        :param SeqRepoDataProxy dp: Access to seqrepo data
        :param UTA uta: Access to UTA queries
        """
        self.dp = dp
        self.uta = uta

    """The Genomic Base class."""
    def get_nc_accessions(self, classification):
        """Get NC accession for a given classification."""
        hgvs = [t.token for t in classification.all_tokens if
                t.token_type in ['HGVS', 'ReferenceSequence']]
        nc_accessions = []
        if hgvs:
            nc_accessions = [hgvs[0].split(':')[0]]
        else:
            gene_tokens = [t for t in classification.all_tokens
                           if t.token_type == 'GeneSymbol']
            if gene_tokens and len(gene_tokens) == 1:
                nc_accessions = self.uta.get_ac_from_gene(gene_tokens[0].token)
        return nc_accessions

    def get_nc_accession(self, identifier):
        """Given an identifier (assembly+chr), return nc accession."""
        nc_accession = None
        try:
            metadata = \
                self.dp.get_metadata(identifier)
        except KeyError:
            logger.warning('Data Proxy unable to get metadata'
                           f'for {identifier}')
        else:
            aliases = [a for a in metadata['aliases'] if
                       a.startswith('refseq:NC_')]
            if aliases:
                nc_accession = aliases[0].split(':')[-1]
        return nc_accession
