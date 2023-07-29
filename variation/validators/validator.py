"""Module for Validation."""
from typing import List, Optional, Tuple
from abc import ABC, abstractmethod

from gene.query import QueryHandler as GeneQueryHandler
from cool_seq_tool.data_sources import SeqRepoAccess, TranscriptMappings, UTADatabase
from cool_seq_tool.schemas import ResidueMode

from variation.schemas.classification_response_schema import (
    Classification,
    ClassificationType,
)
from variation.schemas.token_response_schema import GeneToken
from variation.schemas.validation_response_schema import ValidationResult
from variation.validators.genomic_base import GenomicBase


class Validator(ABC):
    """The validator class."""

    def __init__(
        self,
        seqrepo_access: SeqRepoAccess,
        transcript_mappings: TranscriptMappings,
        uta: UTADatabase,
        gene_normalizer: GeneQueryHandler,
    ) -> None:
        """Initialize the DelIns validator.

        :param seqrepo_access: Access to SeqRepo data
        :param transcript_mappings: Access to transcript mappings
        :param uta: Access to UTA queries
        :param gene_normalizer: Access to gene-normalizer
        """
        self.transcript_mappings = transcript_mappings
        self.seqrepo_access = seqrepo_access
        self.uta = uta
        self.genomic_base = GenomicBase(self.seqrepo_access, self.uta)
        self.gene_normalizer = gene_normalizer

    @abstractmethod
    async def get_accessions(
        self, classification: Classification, errors: List
    ) -> List[str]:
        """Get accessions for a given classification.
        If `classification.nomenclature == Nomenclature.HGVS`, will return the accession
        in the HGVS expression.
        Else, will get all accessions associated to the gene

        :param classification: The classification for list of tokens
        :param errors: List of errors
        :return: List of accessions
        """

    @abstractmethod
    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Check that classification type can be validated by validator.

        :param ClassificationType classification_type: The type of variation
        :return: `True` if classification_type matches validator's
            classification type. `False` otherwise.
        """

    @abstractmethod
    async def get_valid_invalid_results(
        self, classification: Classification, accessions: List
    ) -> List[ValidationResult]:
        """Get list of validation results for a given classification and accessions

        :param classification: A classification for a list of tokens
        :param accessions: A list of accessions for a classification
        :return: List of validation results containing invalid and valid results
        """

    async def validate(self, classification: Classification) -> ValidationResult:
        """Get list of associated accessions for a classification. Use these accessions
        to perform validation checks (pos exists, accession is valid, reference sequence
        matches expected, etc). Gets list of validation results for a given
        classification

        :param classification: A classification for a list of tokens
        :return: List of validation results containing invalid and valid results
        """
        errors = []

        try:
            # NC_ queries do not have gene tokens
            accessions = await self.get_accessions(classification, errors)
        except IndexError:
            accessions = []

        if errors:
            return [
                ValidationResult(
                    accession=None,
                    classification=classification,
                    is_valid=False,
                    errors=errors,
                )
            ]
        validation_results = await self.get_valid_invalid_results(
            classification, accessions
        )
        return validation_results

    def get_protein_accessions(self, gene_token: GeneToken, errors: List) -> List[str]:
        """Get accessions for variations with protein reference sequence.

        :param gene_token: Gene token for a classification
        :param errors: List of errors
        :return: List of possible protein accessions for the variation
        """
        accessions = self.transcript_mappings.protein_transcripts(gene_token.token)
        if not accessions:
            errors.append(
                f"No protein accessions found for gene symbol: {gene_token.token}"
            )
        return accessions

    def get_cdna_accessions(self, gene_token: GeneToken, errors: List) -> List[str]:
        """Get accessions for variations with cDNA reference sequence.

        :param gene_token: Gene token for a classification
        :param errors: List of errors
        :return: List of possible cDNA accessions for the variation
        """
        accessions = self.transcript_mappings.coding_dna_transcripts(gene_token.token)
        if not accessions:
            errors.append(
                f"No cDNA accessions found for gene symbol: {gene_token.token}"
            )
        return accessions

    async def get_genomic_accessions(
        self, classification: Classification, errors: List
    ) -> List[str]:
        """Get genomic RefSeq accessions for variations with genomic reference sequence.

        :param classification: Classification for a list of tokens
        :param errors: List of errors
        :return: List of possible genomic RefSeq accessions for the variation
        """
        accessions = await self.genomic_base.get_nc_accessions(classification)
        if not accessions:
            errors.append("No genomic accessions found")
        return accessions

    async def _validate_gene_pos(
        self,
        gene: str,
        alt_ac: str,
        pos0: int,
        pos1: Optional[int],
        pos2: Optional[int] = None,
        pos3: Optional[int] = None,
        residue_mode: ResidueMode = ResidueMode.RESIDUE,
    ) -> Optional[str]:
        """Validate whether free text genomic query is valid input.
        If invalid input, add error to list of errors

        :param gene: Gene symbol
        :param alt_ac: Genomic accession
        :param pos0: Queried genomic position
        :param pos1: Queried genomic position
        :param pos2: Queried genomic position
        :param pos3: Queried genomic position
        :param residue_mode: Residue mode for positions
        :return: Invalid error message if invalid. Else, `None`
        """
        gene_start_end = {"start": None, "end": None}
        resp = self.gene_normalizer.search(gene, incl="Ensembl")
        if resp.source_matches:
            ensembl_resp = resp.source_matches[0]
            if ensembl_resp.records[0].locations:
                ensembl_loc = ensembl_resp.records[0].locations[0]
                gene_start_end["start"] = ensembl_loc.interval.start.value
                gene_start_end["end"] = ensembl_loc.interval.end.value - 1

        if gene_start_end["start"] is None and gene_start_end["end"] is None:
            return f"gene-normalizer unable to find Ensembl location for gene: {gene}"
        else:
            assembly = await self.uta.get_chr_assembly(alt_ac)
            if assembly:
                # Not in GRCh38 assembly. Gene normalizer only uses 38, so we
                # need to liftover to GRCh37 coords
                chromosome, assembly = assembly
                for key in gene_start_end.keys():
                    gene_pos = gene_start_end[key]
                    gene_pos_liftover = self.uta.liftover_38_to_37.convert_coordinate(
                        chromosome, gene_pos
                    )
                    if gene_pos_liftover is None or len(gene_pos_liftover) == 0:
                        return f"{gene_pos} does not exist on {chromosome}"
                    else:
                        gene_start_end[key] = gene_pos_liftover[0][1]

            gene_start = gene_start_end["start"]
            gene_end = gene_start_end["end"]

            for pos in [pos0, pos1, pos2, pos3]:
                if pos not in ["?", None]:
                    if residue_mode == "residue":
                        pos -= 1
                    if not (gene_start <= pos <= gene_end):
                        return f"Position {pos} out of index on {alt_ac} on gene, {gene}"  # noqa: E501

    def validate_reference_sequence(
        self,
        ac: str,
        start_pos: int,
        end_pos: int,
        expected_ref: str,
        residue_mode: ResidueMode = ResidueMode.RESIDUE,
    ) -> Optional[str]:
        """Validate that expected reference sequence matches actual reference sequence.
        This is also in translator, but there is a ticket to have this method be moved
        to cool-seq-tool. Once added, will be removed

        :param ac: Accession
        :param start_pos: Start position
        :param end_pos: End position
        :param expected_ref: The expected reference sequence (from input string)
        :param residue_mode: Residue mode for `start_pos` and `end_pos`
        :return: Invalid message if invalid. If valid, `None`
        """
        actual_ref, err_msg = self.seqrepo_access.get_reference_sequence(
            ac, start=start_pos, end=end_pos, residue_mode=residue_mode
        )

        if not err_msg and (actual_ref != expected_ref):
            err_msg = (
                f"Expected to find {expected_ref} at positions ({start_pos}, "
                f"{end_pos}) on {ac} but found {actual_ref}"
            )

        return err_msg

    async def get_cds_start(self, ac: str) -> Tuple[Optional[int], Optional[str]]:
        """Get coding start site for accession

        :param ac: Accession to get coding start site for
        :return: Tuple containing coding start site (if successful) and errors
            (if unsuccessful)
        """
        cds_start_end = await self.uta.get_cds_start_end(ac)

        if cds_start_end:
            cds_start = cds_start_end[0]
            msg = None
        else:
            cds_start = None
            msg = f"Unable to get CDS start for accession: {ac}"

        return cds_start, msg

    def validate_ac_and_pos(
        self,
        ac: str,
        start_pos: int,
        end_pos: Optional[int] = None,
        residue_mode: ResidueMode = ResidueMode.RESIDUE,
    ) -> Optional[str]:
        """Validate that accession exists and that position(s) exist on accession

        :param ac: Accession
        :param start_pos: Start position on accession
        :param end_position: End position on accession
        :param residue_mode: Residue mode for `start_pos` and `end_pos`
        :return: If valid accession and position(s) on accession, `None`. If invalid
            accession or invalid position(s) on accession, return error message
        """
        if residue_mode == ResidueMode.RESIDUE:
            start_pos -= 1

        msg = None
        ref_len = None
        try:
            if end_pos:
                ref_len = len(self.seqrepo_access.sr[ac][start_pos:end_pos])
            else:
                ref_len = len(self.seqrepo_access.sr[ac][start_pos])
        except KeyError:
            msg = f"Accession does not exist in SeqRepo: {ac}"
        except ValueError as e:
            msg = f"{e} on accession ({ac})"
        else:
            if end_pos:
                if not ref_len or (end_pos - start_pos != ref_len):
                    msg = f"Positions ({start_pos}, {end_pos}) not valid on accession ({ac})"  # noqa: E501
            else:
                if not ref_len:
                    msg = f"Position ({start_pos}) not valid on accession ({ac})"

        return msg
