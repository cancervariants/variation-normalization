"""Module for Validation."""
from typing import List, Optional, Tuple
from abc import ABC, abstractmethod

from gene.query import QueryHandler as GeneQueryHandler
from cool_seq_tool.data_sources import (
    SeqRepoAccess, TranscriptMappings, MANETranscript, UTADatabase
)
from cool_seq_tool.schemas import ResidueMode

from variation.schemas.classification_response_schema import (
    Classification, ClassificationType
)
from variation.schemas.token_response_schema import GeneToken
from variation.schemas.validation_response_schema import ValidationResult
from variation.tokenizers import GeneSymbol
from variation.validators.genomic_base import GenomicBase
from variation.vrs_representation import VRSRepresentation


class Validator(ABC):
    """The validator class."""

    def __init__(self, seqrepo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol,
                 mane_transcript: MANETranscript,
                 uta: UTADatabase,
                 gene_normalizer: GeneQueryHandler,
                 vrs: VRSRepresentation) -> None:
        """Initialize the DelIns validator.

        :param seqrepo_access: Access to SeqRepo data
        :param transcript_mappings: Access to transcript mappings
        :param gene_symbol: Gene symbol tokenizer
        :param mane_transcript: Access MANE Transcript information
        :param uta: Access to UTA queries
        :param gene_normalizer: Access to gene-normalizer
        :param vrs: Class for creating VRS objects
        """
        self.transcript_mappings = transcript_mappings
        self.seqrepo_access = seqrepo_access
        self._gene_matcher = gene_symbol
        self.uta = uta
        self.genomic_base = GenomicBase(self.seqrepo_access, self.uta)
        self.mane_transcript = mane_transcript
        self.gene_normalizer = gene_normalizer
        self.vrs = vrs

    @abstractmethod
    def variation_name(self) -> str:
        """Return the variation name.

        :return: variation class name
        """

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
            self, classification_type: ClassificationType) -> bool:
        """Check that classification type can be validated by validator.

        :param ClassificationType classification_type: The type of variation
        :return: `True` if classification_type matches validator's
            classification type. `False` otherwise.
        """

    @abstractmethod
    async def get_valid_invalid_results(
        self, classification: Classification, accessions: List
    ) -> None:
        """Add validation result objects to a list of results.

        :param classification: A classification for a list of tokens
        :param accessions: A list of accessions for a classification
        """

    async def validate(
        self, classification: Classification
    ) -> ValidationResult:
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
                    errors=errors
                )
            ]
        validation_results = await self.get_valid_invalid_results(
            classification, accessions
        )
        return validation_results

    def get_protein_accessions(
        self, gene_token: GeneToken, errors: List
    ) -> List[str]:
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

    def get_cdna_accessions(
        self, gene_token: GeneToken, errors: List
    ) -> List[str]:
        """Get accessions for variations with coding DNA reference sequence.

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

    async def _validate_gene_pos(self, gene: str, alt_ac: str, pos1: int, pos2: int,
                                 errors: List, pos3: int = None, pos4: int = None,
                                 residue_mode: str = "residue") -> None:
        """Validate whether free text genomic query is valid input.
        If invalid input, add error to list of errors

        :param str gene: Queried gene
        :param str alt_ac: Genomic accession
        :param int pos1: Queried genomic position
        :param int pos2: Queried genomic position
        :param int pos3: Queried genomic position
        :param int pos4: Queried genomic position
        :param str residue_mode: Must be either `inter-residue` or `residue`
        :param List errors: List of errors
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
            errors.append(f"gene-normalizer unable to find Ensembl location"
                          f"for {gene}")
        else:
            assembly = await self.uta.get_chr_assembly(alt_ac)
            if assembly:
                # Not in GRCh38 assembly. Gene normalizer only uses 38, so we
                # need to liftover to GRCh37 coords
                chromosome, assembly = assembly
                for key in gene_start_end.keys():
                    gene_pos = gene_start_end[key]
                    gene_pos_liftover = self.uta.liftover_38_to_37.convert_coordinate(  # noqa: E501
                        chromosome, gene_pos)
                    if gene_pos_liftover is None or len(gene_pos_liftover) == 0:
                        errors.append(f"{gene_pos} does not"
                                      f" exist on {chromosome}")
                        return None
                    else:
                        gene_start_end[key] = gene_pos_liftover[0][1]

            gene_start = gene_start_end["start"]
            gene_end = gene_start_end["end"]

            for pos in [pos1, pos2, pos3, pos4]:
                if pos not in ["?", None]:
                    if residue_mode == "residue":
                        pos -= 1
                    if not (gene_start <= pos <= gene_end):
                        errors.append(f"Position {pos} out of index on "
                                      f"{alt_ac} on gene, {gene}")

    def _check_index(self, ac: str, pos: int, errors: List) -> Optional[str]:
        """Check that index actually exists

        :param str ac: Accession
        :param int pos: Position changes
        :param List errors: List of errors
        :return: Reference sequence
        """
        seq, w = self.seqrepo_access.get_reference_sequence(ac, pos)
        if not seq:
            errors.append(w)
        return seq

    def validate_reference_sequence(
        self, ac: str, start_pos: int, end_pos: int,
        expected_ref: str, residue_mode: ResidueMode = ResidueMode.RESIDUE
    ) -> Optional[str]:
        """Validate that expected reference sequence matches actual"""
        actual_ref, _ = self.seqrepo_access.get_reference_sequence(
            ac, start=start_pos, end=end_pos, residue_mode=residue_mode
        )

        msg = None
        if actual_ref != expected_ref:
            msg = (f"Expected to find {expected_ref} at positions "
                   f"({start_pos}, {end_pos}) on {ac} but found {actual_ref}")

        return msg

    async def get_cds_start(self, ac: str) -> Tuple[Optional[int], Optional[str]]:
        cds_start_end = await self.uta.get_cds_start_end(ac)

        if cds_start_end:
            cds_start = cds_start_end[0]
            msg = None
        else:
            cds_start = None
            msg = f"Unable to get CDS start for accession: {ac}"

        return cds_start, msg

    def validate_accession(self, ac: str) -> Optional[str]:
        try:
            self.seqrepo_access.sr[ac][0]
        except KeyError:
            msg = f"Accession does not exist in SeqRepo: {ac}"
        else:
            msg = None

        return msg

    def validate_ac_and_pos(
        self, ac: str, start_pos: int, end_pos: Optional[int] = None,
        residue_mode: ResidueMode = ResidueMode.RESIDUE
    ) -> Optional[str]:
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
