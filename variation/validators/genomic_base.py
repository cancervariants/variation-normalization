"""Module for Genomic Validation methods."""
import logging
from typing import List, Optional, Union

from cool_seq_tool.data_sources import SeqRepoAccess, UTADatabase

from variation.schemas.classification_response_schema import (
    Classification,
    Nomenclature,
)
from variation.schemas.service_schema import ClinVarAssembly

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
        self,
        classification: Classification,
        input_assembly: Optional[
            Union[ClinVarAssembly.GRCH37, ClinVarAssembly.GRCH38]
        ] = None,
    ) -> List[str]:
        """Get genomic RefSeq accession for a given classification.

        :param classification: A classification for a list of tokens
        :param input_assembly: Assembly used for initial input query. Only used when
            initial query is using genomic free text or gnomad vcf format
        :return: List of genomic RefSeq accessions
        """
        if classification.nomenclature == Nomenclature.HGVS:
            nc_accessions = [classification.ac]
        elif classification.nomenclature == Nomenclature.FREE_TEXT:
            nc_accessions = await self.uta.get_ac_from_gene(
                classification.gene_token.matched_value
            )

            if input_assembly:
                updated_nc_accessions = []
                for alt_ac in nc_accessions:
                    aliases, _ = self.seqrepo_access.translate_identifier(
                        alt_ac, input_assembly
                    )
                    if aliases:
                        updated_nc_accessions.append(alt_ac)
                        break

                nc_accessions = updated_nc_accessions
        elif classification.nomenclature == Nomenclature.GNOMAD_VCF:
            gnomad_vcf_token = classification.matching_tokens[0]
            chromosome = gnomad_vcf_token.chromosome
            nc_accessions = []

            if input_assembly:
                ac = self.get_nc_accession(f"{input_assembly.value}:{chromosome}")
                if ac:
                    nc_accessions.append(ac)
            else:
                for assembly in [ClinVarAssembly.GRCH38, ClinVarAssembly.GRCH37]:
                    ac = self.get_nc_accession(f"{assembly.value}:{chromosome}")
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
            logger.warning("Data Proxy unable to get metadata" f"for {identifier}")
        else:
            aliases = [a for a in translated_identifiers if a.startswith("refseq:NC_")]
            if aliases:
                nc_accession = aliases[0].split(":")[-1]
        return nc_accession
