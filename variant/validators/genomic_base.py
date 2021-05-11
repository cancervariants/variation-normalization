"""Module for Genomic Validation methods."""
from variant import GENE_NORMALIZER
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
import logging


logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)


class GenomicBase:
    """Genomic Base class for validation methods."""

    def __init__(self, dp: SeqRepoDataProxy):
        """Initialize the Genomic base class."""
        self.dp = dp

    """The Genomic Base class."""
    def get_nc_accessions(self, classification):
        """Get NC accession for a given classification."""
        hgvs = [t.token for t in classification.all_tokens if
                t.token_type in ['HGVS', 'ReferenceSequence']]
        nc_accessions = []
        if hgvs:
            nc_accession = hgvs[0].split(':')[0]
            nc_accessions = \
                self.get_nc_accessions_from_nc_accession(nc_accession)
        else:
            gene_tokens = [t for t in classification.all_tokens
                           if t.token_type == 'GeneSymbol']
            if gene_tokens and len(gene_tokens) == 1:
                resp = GENE_NORMALIZER.search_sources(gene_tokens[0].token,
                                                      incl='hgnc')
                if resp['source_matches'][0]['records']:
                    record = resp['source_matches'][0]['records'][0]
                    loc = record.locations[0] if record.locations else None
                    # TODO: what about multiple chr locations?
                    if loc and loc.chr:
                        for identifier in ['GRCh38', 'GRCh37']:
                            nc_accession =  \
                                self.get_nc_accession(f"{identifier}:"
                                                      f"{loc.chr}")
                            if nc_accession:
                                nc_accessions.append(nc_accession)
        return list(set(nc_accessions))

    def get_nc_accessions_from_nc_accession(self, nc_accession):
        """Given NC accession, find other version from other assembly."""
        nc_accessions = [nc_accession]
        try:
            assembly = None
            for a in self.dp.get_metadata(nc_accession)['aliases']:
                if a.startswith('GRCh3'):
                    assembly = a
                    break
        except KeyError:
            pass
        else:
            if assembly:
                if assembly.startswith('GRCh38'):
                    nc_accession = \
                        self.get_nc_accession(f"GRCh37:"
                                              f"{assembly.split(':')[1]}")
                    if nc_accession:
                        nc_accessions.append(nc_accession)
                elif assembly.startswith('GRCh37'):
                    nc_accession = self.get_nc_accession(
                        f"GRCh38:{assembly.split(':')[1]}")
                    if nc_accession:
                        nc_accessions.append(nc_accession)
        return nc_accessions

    def get_nc_accession(self, identifier):
        """Given an identifier (assembly+chr), return nc accession."""
        nc_accession = None
        try:
            metadata = \
                self.dp.get_metadata(identifier)
        except KeyError:
            logger.warning('Data Proxy unable to get metadata'
                           f'for GRCh38:{identifier}')
        else:
            aliases = [a for a in metadata['aliases'] if
                       a.startswith('refseq:NC_')]
            if aliases:
                nc_accession = aliases[0].split(':')[-1]
        return nc_accession
