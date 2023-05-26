"""Module for Genomic Validation methods."""
import logging
from typing import List, Optional

from cool_seq_tool.data_sources import UTADatabase, SeqRepoAccess

from variation.schemas.classification_response_schema import Classification, Nomenclature
from variation.schemas.token_response_schema import GeneToken, TokenType


logger = logging.getLogger("variation")
logger.setLevel(logging.DEBUG)


class GenomicBase:
    """Genomic Base class for validation methods."""

    def __init__(self, seqrepo_access: SeqRepoAccess, uta: UTADatabase) -> None:
        """Initialize the Genomic base class.

        :param SeqRepoAccess seqrepo_access: Access to seqrepo
        :param UTADatabase uta: Access to UTA queries
        """
        self.seqrepo_access = seqrepo_access
        self.uta = uta

    """The Genomic Base class."""
    async def get_nc_accessions(
        self, classification: Classification, gene_tokens: List[GeneToken]
    ) -> List[str]:
        """Get NC accession for a given classification."""
        hgvs = [t.accession for t in classification.matching_tokens if
                t.token_type == TokenType.HGVS]
        nc_accessions = []
        if hgvs:
            nc_accessions = [hgvs[0]]
        else:
            if classification.nomenclature == Nomenclature.FREE_TEXT:
                if gene_tokens and len(gene_tokens) == 1:
                    nc_accessions = await self.uta.get_ac_from_gene(
                        gene_tokens[0].matched_value
                    )
            else:
                raise NotImplementedError

            # chromosome = [t.token for t in classification.matching_tokens if
            #               t.token_type in [TokenType.CHROMOSOME]]
            # if chromosome:
            #     chromosome = chromosome[0]
            #     for assesmbly in ["GRCh37", "GRCh38"]:
            #         ac = self.get_nc_accession(f"{assesmbly}:{chromosome}")
            #         if ac:
            #             nc_accessions.append(ac)
            # else:
            #     gene_tokens = [t.token for t in classification.matching_tokens
            #                    if t.token_type == TokenType.GENE]
            #     if gene_tokens and len(gene_tokens) == 1:
            #         nc_accessions = await self.uta.get_ac_from_gene(gene_tokens[0])
        return nc_accessions

    def get_nc_accession(self, identifier: str) -> Optional[str]:
        """Given an identifier (assembly+chr), return nc accession."""
        nc_accession = None
        try:
            translated_identifiers, _ = self.seqrepo_access.translate_identifier(
                identifier)
        except KeyError:
            logger.warning("Data Proxy unable to get metadata"
                           f"for {identifier}")
        else:
            aliases = [a for a in translated_identifiers if a.startswith("refseq:NC_")]
            if aliases:
                nc_accession = aliases[0].split(":")[-1]
        return nc_accession
