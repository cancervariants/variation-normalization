"""Module for Validation."""

from abc import ABC, abstractmethod
from typing import Literal

from cool_seq_tool.handlers import SeqRepoAccess
from cool_seq_tool.mappers import LiftOver
from cool_seq_tool.schemas import Assembly, CoordinateType
from cool_seq_tool.sources import TranscriptMappings, UtaDatabase
from gene.query import QueryHandler as GeneQueryHandler
from gene.schemas import SourceName

from variation.schemas.classification_response_schema import (
    AmbiguousType,
    Classification,
    ClassificationType,
    GenomicDeletionAmbiguousClassification,
    GenomicDuplicationAmbiguousClassification,
    ProteinDeletionClassification,
    ProteinDelInsClassification,
    ProteinInsertionClassification,
    ProteinReferenceAgreeClassification,
    ProteinStopGainClassification,
    ProteinSubstitutionClassification,
)
from variation.schemas.token_response_schema import GeneToken
from variation.schemas.validation_response_schema import ValidationResult
from variation.utils import get_aa1_codes
from variation.validators.genomic_base import GenomicBase


class Validator(ABC):
    """The validator class."""

    def __init__(
        self,
        seqrepo_access: SeqRepoAccess,
        transcript_mappings: TranscriptMappings,
        uta: UtaDatabase,
        gene_normalizer: GeneQueryHandler,
        liftover: LiftOver,
    ) -> None:
        """Initialize the DelIns validator.

        :param seqrepo_access: Access to SeqRepo data
        :param transcript_mappings: Access to transcript mappings
        :param uta: Access to UTA queries
        :param gene_normalizer: Access to gene-normalizer
        :param liftover: Instance to provide mapping between human genome assemblies
        """
        self.transcript_mappings = transcript_mappings
        self.seqrepo_access = seqrepo_access
        self.uta = uta
        self.genomic_base = GenomicBase(self.seqrepo_access, self.uta)
        self.gene_normalizer = gene_normalizer
        self.liftover = liftover

    @abstractmethod
    async def get_accessions(
        self, classification: Classification, errors: list
    ) -> list[str]:
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
        self, classification: Classification, accessions: list
    ) -> list[ValidationResult]:
        """Get list of validation results for a given classification and accessions

        :param classification: A classification for a list of tokens
        :param accessions: A list of accessions for a classification
        :return: List of validation results containing invalid and valid results
        """

    async def validate(self, classification: Classification) -> list[ValidationResult]:
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
        return await self.get_valid_invalid_results(classification, accessions)

    def get_protein_accessions(self, gene_token: GeneToken, errors: list) -> list[str]:
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

    def get_cdna_accessions(self, gene_token: GeneToken, errors: list) -> list[str]:
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
        self, classification: Classification, errors: list
    ) -> list[str]:
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
        pos1: int | None,
        pos2: int | None = None,
        pos3: int | None = None,
        coordinate_type: CoordinateType = CoordinateType.RESIDUE,
    ) -> str | None:
        """Validate whether free text genomic query is valid input.
        If invalid input, add error to list of errors

        :param gene: Gene symbol
        :param alt_ac: Genomic accession
        :param pos0: Queried genomic position
        :param pos1: Queried genomic position
        :param pos2: Queried genomic position
        :param pos3: Queried genomic position
        :param coordinate_type: Coordinate type for positions
        :return: Invalid error message if invalid. Else, `None`
        """
        gene_start_end = {"start": None, "end": None}
        resp = self.gene_normalizer.search(gene, incl=SourceName.ENSEMBL.value)
        if resp.source_matches:
            ensembl_resp = resp.source_matches[SourceName.ENSEMBL]
            if all(
                (ensembl_resp, ensembl_resp.records, ensembl_resp.records[0].locations)
            ):
                ensembl_loc = ensembl_resp.records[0].locations[0]
                gene_start_end["start"] = ensembl_loc.start
                gene_start_end["end"] = ensembl_loc.end - 1

        if gene_start_end["start"] is None and gene_start_end["end"] is None:
            return f"gene-normalizer unable to find Ensembl location for gene: {gene}"

        assembly = await self.uta.get_chr_assembly(alt_ac)
        if assembly:
            # Not in GRCh38 assembly. Gene normalizer only uses 38, so we
            # need to liftover to GRCh37 coords
            chromosome, assembly = assembly
            for key in gene_start_end:
                gene_pos = gene_start_end[key]
                gene_pos_liftover = self.liftover.get_liftover(
                    chromosome, gene_pos, Assembly.GRCH37
                )
                if gene_pos_liftover is None or len(gene_pos_liftover) == 0:
                    return f"{gene_pos} does not exist on {chromosome}"
                gene_start_end[key] = gene_pos_liftover[1]

        gene_start = gene_start_end["start"]
        gene_end = gene_start_end["end"]

        for pos in [pos0, pos1, pos2, pos3]:
            if pos not in ["?", None]:
                if coordinate_type == CoordinateType.RESIDUE:
                    pos -= 1
                if not (gene_start <= pos <= gene_end):
                    return f"Position {pos} out of index on {alt_ac} on gene, {gene}"

        return None

    def validate_reference_sequence(
        self,
        ac: str,
        start_pos: int,
        end_pos: int,
        expected_ref: str,
        coordinate_type: CoordinateType = CoordinateType.RESIDUE,
    ) -> str | None:
        """Validate that expected reference sequence matches actual reference sequence.
        This is also in translator, but there is a ticket to have this method be moved
        to cool-seq-tool. Once added, will be removed

        :param ac: Accession
        :param start_pos: Start position
        :param end_pos: End position
        :param expected_ref: The expected reference sequence (from input string)
        :param coordinate_type: Coordinate type for `start_pos` and `end_pos`
        :return: Invalid message if invalid. If valid, `None`
        """
        actual_ref, err_msg = self.seqrepo_access.get_reference_sequence(
            ac, start=start_pos, end=end_pos, coordinate_type=coordinate_type
        )

        if not err_msg and (actual_ref != expected_ref):
            err_msg = (
                f"Expected to find {expected_ref} at positions ({start_pos}, "
                f"{end_pos}) on {ac} but found {actual_ref}"
            )

        return err_msg

    async def get_cds_start(self, ac: str) -> tuple[int | None, str | None]:
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
        end_pos: int | None = None,
        coordinate_type: CoordinateType = CoordinateType.RESIDUE,
    ) -> str | None:
        """Validate that accession exists and that position(s) exist on accession

        :param ac: Accession
        :param start_pos: Start position on accession
        :param end_position: End position on accession
        :param coordinate_type: Coordinate type for `start_pos` and `end_pos`
        :return: If valid accession and position(s) on accession, `None`. If invalid
            accession or invalid position(s) on accession, return error message
        """
        if coordinate_type == CoordinateType.RESIDUE:
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
                    msg = f"Positions ({start_pos}, {end_pos}) not valid on accession ({ac})"
            else:
                if not ref_len:
                    msg = f"Position ({start_pos}) not valid on accession ({ac})"

        return msg

    @staticmethod
    def validate_5_prime_to_3_prime(
        pos0: int,
        pos1: int | Literal["?"] | None,
        pos2: int | Literal["?"] | None = None,
        pos3: int | Literal["?"] | None = None,
    ) -> str | None:
        """Validate that positions are unique and listed from 5' to 3'

        :param pos0: Position 0
        :param pos1: Position 1
        :param pos2: Position 2
        :param pos3: Position 3
        :return: Message if positions are not unique or not listed from 5' to 3'.
            Else, `None`
        """
        prev_pos = None
        invalid_msg = None
        for pos in [pos0, pos1, pos2, pos3]:
            if pos not in {"?", None}:
                if prev_pos is None:
                    prev_pos = pos
                else:
                    if pos <= prev_pos:
                        invalid_msg = (
                            "Positions should contain two different positions and "
                            "should be listed from 5' to 3'"
                        )
                        break

                    prev_pos = pos
        return invalid_msg

    def validate_ambiguous_classification(
        self,
        classification: GenomicDeletionAmbiguousClassification
        | GenomicDuplicationAmbiguousClassification,
    ) -> str | None:
        """Validate that ambiguous type is supported and that positions are unique and
        listed from 5' to 3'

        :param classification: Ambiguous duplication or deletion classification
        :return: Message if ambiguous type is not supported, positions are not unique,
            or if positions are not listed from 5' to 3'. Else, `None`
        """
        invalid_msg = None
        if classification.ambiguous_type not in {
            AmbiguousType.AMBIGUOUS_1,
            AmbiguousType.AMBIGUOUS_2,
            AmbiguousType.AMBIGUOUS_5,
            AmbiguousType.AMBIGUOUS_7,
        }:
            invalid_msg = f"{classification.ambiguous_type} is not yet supported"
        else:
            invalid_msg = self.validate_5_prime_to_3_prime(
                classification.pos0,
                pos1=classification.pos1,
                pos2=classification.pos2,
                pos3=classification.pos3,
            )
        return invalid_msg

    def validate_protein_hgvs_classification(
        self,
        classification: ProteinDelInsClassification
        | ProteinDeletionClassification
        | ProteinInsertionClassification
        | ProteinReferenceAgreeClassification
        | ProteinStopGainClassification
        | ProteinSubstitutionClassification,
    ) -> list[str]:
        """Validate protein HGVS classification

        :param classification: Classification
            Will be mutated if used 3 letter AA codes
        :return: List of invalid error messages if found
        """
        errors = []

        if hasattr(classification, "ref"):
            aa1_ref = get_aa1_codes(classification.ref)
            if aa1_ref:
                classification.ref = aa1_ref
            else:
                errors.append(f"`ref` not valid amino acid(s): {classification.ref}")

        if hasattr(classification, "alt"):
            aa1_alt = get_aa1_codes(classification.alt)
            if aa1_alt:
                classification.alt = aa1_alt
            else:
                errors.append(f"`alt` not valid amino acid(s): {classification.alt}")

        if hasattr(classification, "aa0"):
            aa0_codes = get_aa1_codes(classification.aa0)
            if aa0_codes:
                classification.aa0 = aa0_codes
            else:
                errors.append(f"`aa0` not valid amino acid(s): {classification.aa0}")

        if hasattr(classification, "aa1") and classification.aa1:
            aa1_codes = get_aa1_codes(classification.aa1)
            if aa1_codes:
                classification.aa1 = aa1_codes
            else:
                errors.append(f"`aa1` not valid amino acid(s): {classification.aa1}")

        if hasattr(classification, "inserted_sequence"):
            ins_codes = get_aa1_codes(classification.inserted_sequence)
            if ins_codes:
                classification.inserted_sequence = ins_codes
            else:
                errors.append(
                    f"`inserted_sequence` not valid amino acid(s): "
                    f"{classification.inserted_sequence}"
                )

        return errors
