"""Module for Validation."""
from typing import List, Optional, Dict, Tuple
from abc import ABC, abstractmethod

from ga4gh.vrsatile.pydantic.vrs_models import CopyChange
from gene.query import QueryHandler as GeneQueryHandler
from cool_seq_tool.data_sources import (
    SeqRepoAccess, TranscriptMappings, MANETranscript, UTADatabase
)
from cool_seq_tool.schemas import ResidueMode

from variation.schemas.classification_response_schema import (
    Classification, ClassificationType, Nomenclature
)
from variation.schemas.app_schemas import Endpoint
from variation.schemas.token_response_schema import TokenType, GeneToken
from variation.schemas.validation_response_schema import ValidationResult
from variation.tokenizers import GeneSymbol
from variation.validators.genomic_base import GenomicBase
from variation.schemas.normalize_response_schema import (
    HGVSDupDelMode as HGVSDupDelModeEnum
)
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
        self._protein_ac_to_gene_mapping = [
            self.transcript_mappings.get_gene_symbol_from_ensembl_protein,
            self.transcript_mappings.get_gene_symbol_from_refeq_protein
        ]
        self._cdna_ac_to_gene_mapping = [
            self.transcript_mappings.get_gene_symbol_from_refseq_rna,
            self.transcript_mappings.get_gene_symbol_from_ensembl_transcript
        ]

    @abstractmethod
    def variation_name(self) -> str:
        """Return the variation name.

        :return: variation class name
        """

    @abstractmethod
    def get_gene_tokens(
            self, classification: Classification) -> List[GeneToken]:
        """Return a list of gene tokens for a classification.

        :param Classification classification: Classification for a list of
            tokens
        :return: A list of gene tokens for the classification
        """

    @abstractmethod
    async def get_transcripts(self, gene_tokens: List, classification: Classification,
                              errors: List) -> Optional[List[str]]:
        """Get transcript accessions for a given classification.

        :param List gene_tokens: A list of gene tokens
        :param Classification classification: A classification for a list of
            tokens
        :param List errors: List of errors
        :return: List of transcript accessions
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
        self, classification_tokens: List, transcripts: List,
        classification: Classification, results: List, gene_tokens: List,
        mane_data_found: Dict, is_identifier: bool,
        hgvs_dup_del_mode: HGVSDupDelModeEnum,
        endpoint_name: Optional[Endpoint] = None,
        baseline_copies: Optional[int] = None,
        copy_change: Optional[CopyChange] = None,
        do_liftover: bool = False
    ) -> None:
        """Add validation result objects to a list of results.

        :param List classification_tokens: A list of classification Tokens
        :param List transcripts: A list of transcript accessions
        :param Classification classification: A classification for a list of
            tokens
        :param List results: Stores validation result objects
        :param List gene_tokens: List of GeneMatchTokens for a classification
        :param Dict mane_data_found: MANE Transcript information found
        :param bool is_identifier: `True` if identifier is given for exact
            location. `False` otherwise.
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Must be: `default`,
            `copy_number_count`, `copy_number_change`, `repeated_seq_expr`,
            `literal_seq_expr`. This parameter determines how to represent HGVS dup/del
            expressions as VRS objects.
        :param Optional[Endpoint] endpoint_name: Then name of the endpoint being used
        :param Optional[int] baseline_copies: Baseline copies number
        :param Optional[CopyChange] copy_change: The copy change
        :param bool do_liftover: Whether or not to liftover to GRCh38 assembly
        """

    async def validate(
        self, classification: Classification
    ) -> ValidationResult:
        errors = []

        gene_tokens = self.get_gene_tokens(classification)
        if len(gene_tokens) > 1:
            errors.append("More than one gene symbol found for a single"
                          f" {self.variation_name()}")

        try:
            # NC_ queries do not have gene tokens
            transcripts = await self.get_transcripts(
                gene_tokens, classification, errors
            )
        except IndexError:
            transcripts = []

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
            classification, transcripts, gene_tokens
        )
        return validation_results

    def get_protein_transcripts(self, gene_tokens: List,
                                errors: List) -> Optional[List[str]]:
        """Get transcripts for variations with protein reference sequence.

        :param List gene_tokens: List of gene tokens for a classification
        :param List errors: List of errors
        :return: List of possible transcript accessions for the variation
        """
        transcripts = self.transcript_mappings.protein_transcripts(gene_tokens[0].token)
        if not transcripts:
            errors.append(
                f"No transcripts found for gene symbol {gene_tokens[0].token}"
            )
        return transcripts

    def get_coding_dna_transcripts(self, gene_tokens: List,
                                   errors: List) -> Optional[List[str]]:
        """Get transcripts for variations with coding DNA reference sequence.

        :param List gene_tokens: List of gene tokens for a classification
        :param List errors: List of errors
        :return: List of possible transcript accessions for the variation
        """
        transcripts = self.transcript_mappings.coding_dna_transcripts(
            gene_tokens[0].token
        )
        if not transcripts:
            errors.append(
                f"No transcripts found for gene symbol {gene_tokens[0].token}"
            )
        return transcripts

    async def get_genomic_transcripts(
        self, classification: Classification, gene_tokens: List[GeneToken], errors: List
    ) -> Optional[List[str]]:
        """Get NC accessions for variations with genomic reference sequence.

        :param Classification classification: Classification for a list of
            tokens
        :param List errors: List of errors
        :return: List of possible NC accessions for the variation
        """
        nc_accessions = await self.genomic_base.get_nc_accessions(
            classification, gene_tokens
        )
        if not nc_accessions:
            errors.append("Could not find NC_ accession for {self.variation_name()}")
        return nc_accessions

    @staticmethod
    def get_gene_symbol_tokens(
            classification: Classification) -> List[Optional[GeneToken]]:
        """Return tokens with GeneSymbol token type from a classification.

        :param Classification classification: Classification of input string
        :return: List of Gene Match Tokens
        """
        return [t for t in classification.matching_tokens
                if t.token_type == TokenType.GENE]

    def _add_gene_symbol_to_tokens(self, gene_symbol: str, gene_symbols: List,
                                   gene_tokens: List) -> None:
        """Add a gene symbol to list of gene match tokens.

        :param str gene_symbol: Gene symbol
        :param List gene_symbols: List of gene symbols matched
        :param List gene_tokens: List of GeneMatchTokens
        """
        if gene_symbol and gene_symbol not in gene_symbols:
            gene_symbols.append(gene_symbol)
            gene_tokens.append(self._gene_matcher.match(
                gene_symbol))

    def _get_gene_tokens(self, classification: Classification,
                         mappings: List) -> List[Optional[GeneToken]]:
        """Get gene symbol tokens for protein or transcript reference
        sequences.

        :param Classification classification: Classification for a list of
            tokens
        :param List mappings: List of transcript mapping methods for
            corresponding reference sequence
        :return: A list of gene match tokens
        """
        gene_tokens = self.get_gene_symbol_tokens(classification)
        if not gene_tokens:
            if classification.nomenclature == Nomenclature.HGVS:
                accession = classification.ac
            else:
                raise NotImplementedError

            if not accession:
                return []

            gene_symbols = list()
            for mapping in mappings:
                gene_symbol = mapping(accession)
                self._add_gene_symbol_to_tokens(
                    gene_symbol, gene_symbols, gene_tokens
                )
                if gene_tokens:
                    break

        return gene_tokens

    def get_protein_gene_symbol_tokens(
        self, classification: Classification
    ) -> List[GeneToken]:
        """Return gene tokens for a classification with protein reference
        sequence.

        :param Classification classification: The classification for a list of
            tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self._get_gene_tokens(classification, self._protein_ac_to_gene_mapping)

    def get_coding_dna_gene_symbol_tokens(
            self, classification: Classification
    ) -> List[GeneToken]:
        """Return gene symbol tokens for classifications with coding dna
        reference sequence.

        :param Classification classification: Classification of input string
        :return: A list of gene match tokens
        """
        return self._get_gene_tokens(classification, self._cdna_ac_to_gene_mapping)

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
    ):
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
                    msg = f"Positions ({start_pos}, {end_pos}) not valid on accession ({ac})"
            else:
                if not ref_len:
                    msg = f"Position ({start_pos}) not valid on accession ({ac})"

        return msg
