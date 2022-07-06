"""Module for Validation."""
import copy
from typing import List, Optional, Dict, Tuple
from abc import ABC, abstractmethod
import logging

from ga4gh.vrsatile.pydantic.vrs_models import RelativeCopyClass
from gene.query import QueryHandler as GeneQueryHandler
from ga4gh.vrs.extras.translator import Translator
from uta_tools.data_sources import SeqRepoAccess, TranscriptMappings, MANETranscript, \
    UTADatabase

from variation.schemas.classification_response_schema import Classification, \
    ClassificationType
from variation.schemas.app_schemas import Endpoint
from variation.schemas.token_response_schema import GeneMatchToken, Token, \
    GenomicSubstitutionToken
from variation.schemas.validation_response_schema import ValidationResult
from variation.tokenizers import GeneSymbol
from variation.validators.genomic_base import GenomicBase
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum
from variation.vrs_representation import VRSRepresentation

logger = logging.getLogger("variation")
logger.setLevel(logging.DEBUG)


class Validator(ABC):
    """The validator class."""

    def __init__(self, seqrepo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol,
                 mane_transcript: MANETranscript,
                 uta: UTADatabase, tlr: Translator,
                 gene_normalizer: GeneQueryHandler,
                 vrs: VRSRepresentation) -> None:
        """Initialize the DelIns validator.

        :param SeqRepoAccess seqrepo_access: Access to SeqRepo data
        :param TranscriptMappings transcript_mappings: Access to transcript
            mappings
        :param GeneSymbol gene_symbol: Gene symbol tokenizer
        :param MANETranscript mane_transcript: Access MANE Transcript
            information
        :param UTADatabase uta: Access to UTA queries
        :param Translator tlr: Class for translating nomenclatures to and from VRS
        :param GeneQueryHandler gene_normalizer: Access to gene-normalizer
        :param VRSRepresentation vrs: Class for creating VRS objects
        """
        self.transcript_mappings = transcript_mappings
        self.seqrepo_access = seqrepo_access
        self._gene_matcher = gene_symbol
        self.tlr = tlr
        self.uta = uta
        self.genomic_base = GenomicBase(self.seqrepo_access, self.uta)
        self.mane_transcript = mane_transcript
        self.gene_normalizer = gene_normalizer
        self.vrs = vrs

    @abstractmethod
    def is_token_instance(self, t: Token) -> bool:
        """Check to see if token is instance of a token type.

        :param Token t: Classification token to find type of
        :return: `True` if token is instance of class token. `False` otherwise.
        """
        raise NotImplementedError

    @abstractmethod
    def variation_name(self) -> str:
        """Return the variation name.

        :return: variation class name
        """
        raise NotImplementedError

    @abstractmethod
    def get_gene_tokens(
            self, classification: Classification) -> List[GeneMatchToken]:
        """Return a list of gene tokens for a classification.

        :param Classification classification: Classification for a list of
            tokens
        :return: A list of gene tokens for the classification
        """
        raise NotImplementedError

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
        raise NotImplementedError

    @abstractmethod
    def validates_classification_type(
            self, classification_type: ClassificationType) -> bool:
        """Check that classification type can be validated by validator.

        :param ClassificationType classification_type: The type of variation
        :return: `True` if classification_type matches validator's
            classification type. `False` otherwise.
        """
        raise NotImplementedError

    @abstractmethod
    async def get_valid_invalid_results(
        self, classification_tokens: List, transcripts: List,
        classification: Classification, results: List, gene_tokens: List,
        mane_data_found: Dict, is_identifier: bool,
        hgvs_dup_del_mode: HGVSDupDelModeEnum,
        endpoint_name: Optional[Endpoint] = None,
        baseline_copies: Optional[int] = None,
        relative_copy_class: Optional[RelativeCopyClass] = None,
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
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Must be: `default`, `absolute_cnv`,
            `relative_cnv`, `repeated_seq_expr`, `literal_seq_expr`. This parameter
            determines how to represent HGVS dup/del expressions as VRS objects.
        :param Optional[Endpoint] endpoint_name: Then name of the endpoint being used
        :param Optional[int] baseline_copies: Baseline copies number
        :param Optional[RelativeCopyClass] relative_copy_class: The relative copy class
        :param bool do_liftover: Whether or not to liftover to GRCh38 assembly
        """
        raise NotImplementedError

    async def validate(
            self, classification: Classification,
            hgvs_dup_del_mode: HGVSDupDelModeEnum = HGVSDupDelModeEnum.DEFAULT,
            endpoint_name: Optional[Endpoint] = None,
            baseline_copies: Optional[int] = None,
            relative_copy_class: Optional[RelativeCopyClass] = None,
            do_liftover: bool = False
    ) -> List[ValidationResult]:
        """Return validation result for a given classification.

        :param Classification classification: A classification for a list of
            tokens
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Must be: `default`, `absolute_cnv`,
            `relative_cnv`, `repeated_seq_expr`, `literal_seq_expr`. This parameter
            determines how to represent HGVS dup/del expressions as VRS objects.
        :param Optional[Endpoint] endpoint_name: Then name of the endpoint being used
        :param Optional[int] baseline_copies: Baseline copies number
        :param Optional[RelativeCopyClass] relative_copy_class: The relative copy class
        :param bool do_liftover: Whether or not to liftover to GRCh38 assembly
        :return: List of ValidationResult's containing valid and invalid
            results
        """
        results = list()
        errors = list()

        classification_tokens = self.get_classification_tokens(classification)
        if len(classification.non_matching_tokens) > 0:
            errors.append(f"Non matching tokens found for "
                          f"{self.variation_name()}.")

        gene_tokens = self.get_gene_tokens(classification)
        if len(gene_tokens) > 1:
            errors.append("More than one gene symbol found for a single"
                          f" {self.variation_name()}")

        try:
            # NC_ queries do not have gene tokens
            transcripts = await self.get_transcripts(
                gene_tokens, classification, errors)
        except IndexError:
            transcripts = list()

        if len(errors) > 0:
            return [
                self.get_validation_result(
                    classification, None, False, 0, {}, errors, gene_tokens
                )
            ]

        mane_data_found = {
            "mane_select": dict(),
            "mane_plus_clinical": dict(),
            "longest_compatible_remaining": dict(),
            "grch38": dict()
        }

        # If is_identifier, should only run once
        if "HGVS" in classification.matching_tokens:
            is_identifier = True
        else:
            is_identifier = False

        await self.get_valid_invalid_results(
            classification_tokens, transcripts, classification,
            results, gene_tokens, mane_data_found,
            is_identifier, hgvs_dup_del_mode, endpoint_name, baseline_copies,
            relative_copy_class, do_liftover
        )
        return results

    @staticmethod
    def get_validation_result(
            classification: Classification, classification_token: Token,
            is_valid: bool, confidence_score: int, variation: Dict,
            errors: List, gene_tokens: List, identifier: str = None,
            is_mane_transcript: bool = False) -> ValidationResult:
        """Return a validation result object.

        :param Classification classification: The classification for tokens
        :param Token classification_token: Classification token
        :param bool is_valid: Whether or not the classification is valid
        :param int confidence_score: The classification confidence score
        :param Dict variation: A VRS Variation object
        :param List errors: A list of errors for the classification
        :param List gene_tokens: List of GeneMatchTokens
        :param str identifier: Identifier for variation
        :param bool is_mane_transcript: `True` if result is MANE transcript.
            `False` otherwise.
        :return: A validation result
        """
        return ValidationResult(
            classification=classification,
            classification_token=classification_token,
            is_valid=is_valid,
            confidence_score=confidence_score,
            variation=variation,
            errors=errors,
            gene_tokens=gene_tokens,
            is_mane_transcript=is_mane_transcript,
            identifier=identifier
        )

    def get_protein_transcripts(self, gene_tokens: List,
                                errors: List) -> Optional[List[str]]:
        """Get transcripts for variations with protein reference sequence.

        :param List gene_tokens: List of gene tokens for a classification
        :param List errors: List of errors
        :return: List of possible transcript accessions for the variation
        """
        transcripts = self.transcript_mappings.protein_transcripts(gene_tokens[0].token)
        if not transcripts:
            errors.append(f"No transcripts found for gene symbol "
                          f"{gene_tokens[0].token}")
        return transcripts

    def get_coding_dna_transcripts(self, gene_tokens: List,
                                   errors: List) -> Optional[List[str]]:
        """Get transcripts for variations with coding DNA reference sequence.

        :param List gene_tokens: List of gene tokens for a classification
        :param List errors: List of errors
        :return: List of possible transcript accessions for the variation
        """
        transcripts = self.transcript_mappings.coding_dna_transcripts(
            gene_tokens[0].token)
        if not transcripts:
            errors.append(f"No transcripts found for gene symbol "
                          f"{gene_tokens[0].token}")
        return transcripts

    async def get_genomic_transcripts(self, classification: Classification,
                                      errors: List) -> Optional[List[str]]:
        """Get NC accessions for variations with genomic reference sequence.

        :param Classification classification: Classification for a list of
            tokens
        :param List errors: List of errors
        :return: List of possible NC accessions for the variation
        """
        nc_accessions = await self.genomic_base.get_nc_accessions(classification)
        if not nc_accessions:
            errors.append("Could not find NC_ accession for "
                          f"{self.variation_name()}")
        return nc_accessions

    def get_classification_tokens(
            self, classification: Classification
    ) -> List[Optional[Classification]]:
        """Get classification tokens for a given instance.

        :param Classification classification: A classification for a list of
            tokens
        :return: A list of classification tokens
        """
        return [t for t in classification.all_tokens
                if self.is_token_instance(t)]

    @staticmethod
    def get_gene_symbol_tokens(
            classification: Classification) -> List[Optional[GeneMatchToken]]:
        """Return tokens with GeneSymbol token type from a classification.

        :param Classification classification: Classification of input string
        :return: List of Gene Match Tokens
        """
        return [t for t in classification.all_tokens
                if t.token_type == "GeneSymbol"]

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
                         mappings: List) -> List[Optional[GeneMatchToken]]:
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
            refseq = \
                ([t.token for t in classification.all_tokens if
                  t.token_type in ["HGVS", "ReferenceSequence",
                                   "LocusReferenceGenomic"]] or [None])[0]

            if not refseq:
                return []

            if ":" in refseq:
                refseq = refseq.split(":")[0]

            gene_symbols = list()
            for mapping in mappings:
                gene_symbol = mapping(refseq)
                self._add_gene_symbol_to_tokens(
                    gene_symbol, gene_symbols, gene_tokens
                )
                if gene_tokens:
                    break
        return gene_tokens

    def get_protein_gene_symbol_tokens(
            self, classification: Classification
    ) -> List[Optional[GeneMatchToken]]:
        """Return gene tokens for a classification with protein reference
        sequence.

        :param Classification classification: The classification for a list of
            tokens
        :return: A list of Gene Match Tokens in the classification
        """
        mappings = [
            self.transcript_mappings.get_gene_symbol_from_ensembl_protein,
            self.transcript_mappings.get_gene_symbol_from_refeq_protein,
            self.transcript_mappings.get_gene_symbol_from_lrg
        ]
        return self._get_gene_tokens(classification, mappings)

    def get_coding_dna_gene_symbol_tokens(
            self, classification: Classification
    ) -> List[Optional[GeneMatchToken]]:
        """Return gene symbol tokens for classifications with coding dna
        reference sequence.

        :param Classification classification: Classification of input string
        :return: A list of gene match tokens
        """
        mappings = [
            self.transcript_mappings.get_gene_symbol_from_refseq_rna,
            self.transcript_mappings.get_gene_symbol_from_ensembl_transcript,  # noqa: E501
            self.transcript_mappings.get_gene_symbol_from_lrg
        ]
        return self._get_gene_tokens(classification, mappings)

    @staticmethod
    def get_accession(t: str, classification: Classification) -> str:
        """Return accession for a classification

        :param str t: Accession
        :param Classification classification: Classification for token
        :return: Accession
        """
        if "HGVS" in classification.matching_tokens or \
                "ReferenceSequence" in classification.matching_tokens:
            hgvs_token = [t for t in classification.all_tokens if
                          isinstance(t, Token) and t.token_type
                          in ["HGVS", "ReferenceSequence"]][0]
            hgvs_expr = hgvs_token.input_string
            t = hgvs_expr.split(":")[0]
        return t

    def add_validation_result(
            self, variation: Dict, valid_variations: List, results: List,
            classification: Classification, s: Token, t: str,
            gene_tokens: List, errors: List, identifier: str = None,
            is_mane_transcript: bool = False) -> bool:
        """Add validation result to list of results.

        :param Dict variation: A VRS Variation object
        :param List valid_variations: A list containing current valid
            variations
        :param List results: A list of validation results
        :param Classification classification: The classification for tokens
        :param Token s: The classification token
        :param string t: Transcript
        :param List gene_tokens: List of GeneMatchTokens
        :param List errors: A list of errors for the classification
        :param str identifier: Identifier for variation
        :param bool is_mane_transcript: `True` if result is MANE transcript.
            `False` otherwise.
        """
        if not errors:
            if is_mane_transcript or \
                    (variation and variation not in valid_variations):
                results.append(
                    self.get_validation_result(
                        classification, s, True, 1, variation, [], gene_tokens,
                        identifier=identifier if identifier else t,
                        is_mane_transcript=is_mane_transcript
                    )
                )
                valid_variations.append(variation)
                return True
        else:
            results.append(
                self.get_validation_result(
                    classification, s, False, 1, variation, errors, gene_tokens,
                    identifier=identifier if identifier else t,
                    is_mane_transcript=is_mane_transcript
                )
            )
            return False

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

    @staticmethod
    def _get_coord_alt(coordinate: str, mane: Dict,
                       s_copy: Token) -> Optional[Tuple[str, str]]:
        """Get coordinate and alteration

        :param str coordinate: Coordinate used. Must be either `p`, `c`, or `g`
        :param Dict mane: Mane data
        :param Token s_copy: classification token
        :return: Coordinate, alteration
        """
        if coordinate == "g" and mane["status"].lower() != "grch38":
            s_copy.molecule_context = "transcript"
            s_copy.coordinate_type = "c"
            coordinate = s_copy.coordinate_type

            if isinstance(s_copy, GenomicSubstitutionToken) and \
                    mane["strand"] == "-":
                ref_rev = s_copy.ref_nucleotide[::-1]
                alt_rev = s_copy.new_nucleotide[::-1]

                complements = {
                    "A": "T",
                    "T": "A",
                    "C": "G",
                    "G": "C"
                }

                s_copy.ref_nucleotide = ""
                s_copy.new_nucleotide = ""
                for nt in ref_rev:
                    s_copy.ref_nucleotide += complements[nt]
                for nt in alt_rev:
                    s_copy.new_nucleotide += complements[nt]
                alt = s_copy.new_nucleotide
            else:
                alt = None
            return coordinate, alt
        return None

    def add_mane_data(
            self, mane: Dict, mane_data: Dict, coordinate: str, alt_type: str,
            s: Token, alt: str = None, mane_variation: Dict = None) -> None:
        """Add mane transcript information to mane_data.

        :param Dict mane: MANE data
        :param Dict mane_data: All MANE data found for given query
        :param str coordinate: Coordinate used. Must be either `p`, `c`, or `g`
        :param str alt_type: Type of alteration
        :param Token s: Classification token
        :param str alt: Alteration
        :param Dict mane_variation: VRS Variation for mane data
        """
        if not mane:
            return None

        s_copy = copy.deepcopy(s)
        coord_alt = self._get_coord_alt(coordinate, mane, s_copy)
        if coord_alt:
            coordinate = coord_alt[0] if coord_alt[0] else coordinate
            alt = coord_alt[1] if coord_alt[1] else alt
        if mane_variation is None:
            new_allele = self.vrs.to_vrs_allele(
                mane["refseq"], mane["pos"][0] + 1, mane["pos"][1] + 1,
                coordinate, alt_type, [],
                cds_start=mane.get("coding_start_site", None), alt=alt
            )
            variation = new_allele
        else:
            variation = mane_variation

        if not variation:
            return None

        self._add_dict_to_mane_data(mane["refseq"], s_copy, variation,
                                    mane_data, mane["status"])

    @staticmethod
    def _add_dict_to_mane_data(ac: str, s: Token, variation: Dict,
                               mane_data: Dict, status: str) -> None:
        """Add variation data to mane data for normalize endpoint.

        :param str ac: Accession
        :param Token s: Classification token
        :param Dict variation: VRS Variation object
        :param Dict mane_data: MANE Transcript data found for given query
        :param str status: Status for variation (GRCh38, MANE Select,
            MANE Clinical Plus)
        """
        _id = variation["_id"]
        key = "_".join(status.lower().split())

        if _id in mane_data[key].keys():
            mane_data[key][_id]["count"] += 1
        else:
            mane_data[key][_id] = {
                "classification_token": s,
                "accession": ac,
                "count": 1,
                "variation": variation,
                "label": ac  # TODO: Use VRS to translate
            }

    def add_mane_to_validation_results(
            self, mane_data: Dict, valid_alleles: List, results: List,
            classification: Classification, gene_tokens: List) -> None:
        """Add MANE Transcript data to list of validation results.

        :param Dict mane_data: MANE Transcript data found for given query
        :param List valid_alleles: A list containing current valid alleles
        :param List results: A list of validation results
        :param Classification classification: The classification for tokens
        :param List gene_tokens: List of GeneMatchTokens
        """
        mane_data_keys = mane_data.keys()
        for key in ["mane_select", "mane_plus_clinical", "grch38",
                    "longest_compatible_remaining"]:
            highest_count = 0
            mane_result = None
            mane_allele = None
            identifier = None
            if key not in mane_data_keys:
                continue
            _mane_data_keys = mane_data[key].keys()
            for mane_allele_id in _mane_data_keys:
                data = mane_data[key][mane_allele_id]

                if data["count"] > highest_count:
                    highest_count = data["count"]
                    mane_result = data
                    mane_allele = data["variation"]
                    identifier = data["accession"]

            if mane_allele:
                self.add_validation_result(
                    mane_allele, valid_alleles, results, classification,
                    mane_result["classification_token"],
                    mane_result["accession"], gene_tokens, [],
                    identifier=identifier, is_mane_transcript=True
                )
                return

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

    def _is_grch38_assembly(self, t: str) -> bool:
        """Return whether or not accession is GRCh38 assembly.

        :param str t: Accession
        :return: `True` if accession is GRCh38 assembly. `False` otherwise
        """
        translated_identifiers, w = self.seqrepo_access.translate_identifier(t)
        if translated_identifiers:
            return "GRCh38" in ([a for a in translated_identifiers if a.startswith("GRCh")] or [None])[0]  # noqa: E501
        return False

    async def add_genomic_liftover_to_results(
        self, grch38: Dict, errors: List, alt: str, valid_alleles: List, results: List,
        classification: Classification, s: Token, t: str, gene_tokens: List
    ) -> None:
        """Add genomic liftover data to results if genomic GRCh38 variation found
        Currently only used for to_canonical_variation endpoint with genomic data

        :param Dict grch38: GRCh38 data containing accession and position
        :param List errors: List of errors
        :param str alt: Altered sequence
        :param List valid_alleles: List of valid alleles
        :param List results: List of results data
        :param Classification classification: A classification for a list of tokens
        :param Token s: Classification token
        :param str t: Accession
        :param List gene_tokens: List of GeneMatchTokens for a classification
        """
        if grch38:
            self._check_index(grch38["ac"], grch38["pos"][0], errors)
            if grch38["pos"][0] != grch38["pos"][1]:
                self._check_index(grch38["ac"], grch38["pos"][1], errors)

            if not errors:
                variation = self.vrs.to_vrs_allele(
                    grch38["ac"], grch38["pos"][0], grch38["pos"][1],
                    s.coordinate_type, s.alt_type, errors, alt=alt)
                if variation:
                    self.add_validation_result(
                        variation, valid_alleles, results, classification, s, t,
                        gene_tokens, errors, identifier=grch38["ac"],
                        is_mane_transcript=True)
