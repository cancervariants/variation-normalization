"""Module for Genomic Validation methods."""
import logging
from typing import List, Optional

from cool_seq_tool.handlers import SeqRepoAccess
from cool_seq_tool.sources import UtaDatabase

from variation.schemas.classification_response_schema import (
    Classification,
    Nomenclature,
)

_logger = logging.getLogger(__name__)


class GenomicBase:
    """Genomic Base class for validation methods."""

    def __init__(self, seqrepo_access: SeqRepoAccess, uta: UtaDatabase) -> None:
        """Initialize the Genomic base class.

        :param SeqRepoAccess seqrepo_access: Access to seqrepo
        :param UtaDatabase uta: Access to UTA queries
        """
        self.seqrepo_access = seqrepo_access
        self.uta = uta

    """The Genomic Base class."""

    async def get_nc_accessions(self, classification: Classification) -> List[str]:
        """Get NC accession for a given classification."""
        if classification.nomenclature == Nomenclature.HGVS:
            nc_accessions = [classification.ac]
        elif classification.nomenclature == Nomenclature.FREE_TEXT:
            nc_accessions = await self.uta.get_ac_from_gene(
                classification.gene_token.matched_value
            )
        elif classification.nomenclature == Nomenclature.GNOMAD_VCF:
            gnomad_vcf_token = classification.matching_tokens[0]
            chromosome = gnomad_vcf_token.chromosome
            nc_accessions = []

            for assembly in ["GRCh37", "GRCh38"]:
                ac = self.get_nc_accession(f"{assembly}:{chromosome}")
                if ac:
                    nc_accessions.append(ac)
        else:
            raise NotImplementedError

        return nc_accessions

    def get_nc_accession(self, identifier: str) -> Optional[str]:
        """Given an identifier (assembly+chr), return nc accession."""
        nc_accession = None
        try:
            translated_identifiers, _ = self.seqrepo_access.translate_identifier(
                identifier
            )
        except KeyError:
            _logger.warning("Data Proxy unable to get metadata for %s", identifier)
        else:
            aliases = [a for a in translated_identifiers if a.startswith("refseq:NC_")]
            if aliases:
                nc_accession = aliases[0].split(":")[-1]
        return nc_accession
