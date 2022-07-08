"""The base class for Duplication and Deletion Validation."""
import logging
from typing import Optional, List, Dict, Tuple

from ga4gh.vrsatile.pydantic.vrs_models import RelativeCopyClass
from ga4gh.vrs.extras.translator import Translator
from gene.query import QueryHandler as GeneQueryHandler
from uta_tools.data_sources import SeqRepoAccess, TranscriptMappings, UTADatabase, \
    MANETranscript

from variation.schemas.classification_response_schema import Classification, \
    ClassificationType
from variation.schemas.app_schemas import Endpoint
from variation.schemas.token_response_schema import Token, GeneMatchToken
from variation.validators.validator import Validator
from variation.hgvs_dup_del_mode import HGVSDupDelMode
from variation.tokenizers import GeneSymbol
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum
from variation.vrs_representation import VRSRepresentation

logger = logging.getLogger("variation")
logger.setLevel(logging.DEBUG)


class DuplicationDeletionBase(Validator):
    """The Deletion Validator Base class."""

    def __init__(self, seq_repo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol,
                 mane_transcript: MANETranscript,
                 uta: UTADatabase, tlr: Translator,
                 gene_normalizer: GeneQueryHandler, vrs: VRSRepresentation) -> None:
        """Initialize the Deletion Base validator.

        :param SeqRepoAccess seq_repo_access: Access to SeqRepo data
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
        super().__init__(
            seq_repo_access, transcript_mappings, gene_symbol, mane_transcript,
            uta, tlr, gene_normalizer, vrs
        )
        self.hgvs_dup_del_mode = HGVSDupDelMode(seq_repo_access)

    def is_token_instance(self, t: Token) -> bool:
        """Check to see if token is instance of a token type.

        :param Token t: Classification token to find type of
        :return: `True` if token is instance of class token. `False` otherwise.
        """
        raise NotImplementedError

    def variation_name(self) -> str:
        """Return the variation name.

        :return: variation class name
        """
        raise NotImplementedError

    def get_gene_tokens(self, classification: Classification) -> List[GeneMatchToken]:
        """Return a list of gene tokens for a classification.

        :param Classification classification: Classification for a list of
            tokens
        :return: A list of gene tokens for the classification
        """
        raise NotImplementedError

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

    def validates_classification_type(
            self, classification_type: ClassificationType) -> bool:
        """Check that classification type can be validated by validator.

        :param ClassificationType classification_type: The type of variation
        :return: `True` if classification_type matches validator's
            classification type. `False` otherwise.
        """
        raise NotImplementedError

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

    def get_reference_sequence(
            self, ac: str, start: int, end: int, errors: List,
            cds_start: int = None) -> Optional[str]:
        """Get deleted reference sequence.

        :param str ac: Accession
        :param int start: Start position
        :param int end: End position
        :param List errors: List of errors
        :param int cds_start: Coding start site
        :return: Reference sequence of nucleotides
        """
        if cds_start:
            start += cds_start
            if end is not None:
                end += cds_start

        if start and not end:
            ref_sequence, _ = self.seqrepo_access.get_reference_sequence(ac, start)
        elif start is not None and end is not None:
            ref_sequence, _ = self.seqrepo_access.get_reference_sequence(ac, start,
                                                                         end + 1)
        else:
            ref_sequence = None

        if not ref_sequence:
            errors.append("Unable to get reference sequence.")
        return ref_sequence

    def check_reference_sequence(self, t: str, s: Token, errors: List,
                                 cds_start: int = None) -> bool:
        """Check that reference sequence matches deleted sequence.

        :param str t: Accession
        :param Token s: Classification token
        :param List errors: List of errors
        :param int cds_start: Coding start site
        :return: `True` if ref_sequences matches deleted_sequence.
            `False` otherwise.
        """
        ref_sequence = self.get_reference_sequence(
            t, s.start_pos_del, s.end_pos_del, errors, cds_start=cds_start)

        if errors:
            return False
        else:
            if ref_sequence and s.deleted_sequence:
                if ref_sequence != s.deleted_sequence:
                    errors.append(f"Expected deleted sequence {ref_sequence} "
                                  f"but got {s.deleted_sequence}")
                    return False
        return True

    async def add_normalized_genomic_dup_del(
            self, s: Token, t: str, start: int, end: int, gene: str,
            so_id: str, errors: List, hgvs_dup_del_mode: HGVSDupDelModeEnum,
            mane_data_found: Dict, baseline_copies: Optional[int] = None,
            relative_copy_class: Optional[RelativeCopyClass] = None) -> None:
        """Add normalized genomic dup or del to mane data

        :param Token s: Classification token
        :param str t: Accession
        :param int start: Start position
        :param int end: ENd position
        :param str gene: Gene
        :param str so_id: Sequence ontology id
        :param List errors: List of errors
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Must be: `default`, `absolute_cnv`,
            `relative_cnv`, `repeated_seq_expr`, `literal_seq_expr`. This parameter
            determines how to represent HGVS dup/del expressions as VRS objects.
        :param Dict mane_data_found: MANE Transcript information found
        :param Optional[int] baseline_copies: Baseline copies number
        :param Optional[RelativeCopyClass] relative_copy_class: The relative copy class
        """
        mane = await self.mane_transcript.get_mane_transcript(
            t, start, s.coordinate_type, end_pos=end, gene=gene,
            try_longest_compatible=True)

        if mane:
            s.coordinate_type = "c"
            s.molecule_context = "transcript"
            s.so_id = so_id

            allele = self.vrs.to_vrs_allele(
                mane["refseq"], mane["pos"][0] + 1, mane["pos"][1] + 1,
                s.coordinate_type, s.alt_type, errors,
                cds_start=mane["coding_start_site"])

            mane_variation = self.hgvs_dup_del_mode.interpret_variation(
                s.alt_type, allele, errors, hgvs_dup_del_mode,
                baseline_copies=baseline_copies,
                relative_copy_class=relative_copy_class)

            if mane_variation:
                self._add_dict_to_mane_data(mane["refseq"], s, mane_variation,
                                            mane_data_found, mane["status"])

    async def validate_gene_or_accession_pos(
        self, t: str, pos_list: List, errors: List, gene: Optional[str] = None
    ) -> None:
        """Validate positions on gene or accession.
        If not valid, add to list of errors

        :param str t: Accession
        :param List pos_list: List of positions to validate
        :param List errors: List of errors
        :param Optional[str] gene: Gene
        """
        if gene:
            len_pos_list = len(pos_list)
            pos1, pos2 = pos_list[0], pos_list[1]
            if len_pos_list == 2:
                pos3, pos4 = None, None
            elif len_pos_list == 4:
                pos3, pos4 = pos_list[2], pos_list[3]
            else:
                errors.append(f"Unexpected amount of positions: {len_pos_list}")
                return
            await self._validate_gene_pos(gene, t, pos1, pos2, errors, pos3=pos3,
                                          pos4=pos4)
        else:
            for pos in pos_list:
                self._check_index(t, pos, errors)

    async def get_grch38_pos_ac(
            self, t: str, pos1: int, pos2: int,
            pos3: Optional[int] = None, pos4: Optional[int] = None
    ) -> Tuple[str, int, int, Optional[int], Optional[int], Optional[Dict]]:
        """Get GRCh38 start/end and accession

        :param str t: Accession
        :param int pos1: Position 1 (assumes 1-based)
        :param int pos2: Position 2 (assumes 1-based)
        :param int pos3: Position 3 (assumes 1-based)
        :param int pos4: Position 4 (assumes 1-based)
        :return: GRCh38 data (accession and positions)
        """
        grch38_start = await self.mane_transcript.g_to_grch38(t, pos1, pos2)
        if grch38_start:
            pos1, pos2 = grch38_start["pos"]
            if pos3 is not None and pos4 is not None:
                grch38_end = await self.mane_transcript.g_to_grch38(t, pos3, pos4)
                if grch38_end:
                    pos3, pos4 = grch38_end["pos"]
            t = grch38_start["ac"]
        return t, pos1, pos2, pos3, pos4, grch38_start

    def add_grch38_to_mane_data(
            self, t: str, s: Token, errors: List,
            grch38: Dict, mane_data_found: Dict,
            hgvs_dup_del_mode: HGVSDupDelModeEnum,
            ival: Optional[Tuple] = None,
            use_vrs_allele_range: bool = True,
            baseline_copies: Optional[int] = None,
            relative_copy_class: Optional[RelativeCopyClass] = None) -> None:
        """Add grch38 variation to mane data

        :param str t: Accession
        :param Token s: Classification token
        :param List errors: List of errors
        :param Dict grch38: GRCh38 data
        :param Dict mane_data_found: MANE data found for initial query
        :param HGVSDupDelModeEnum hgvs_dup_del_mode: Must be: `default`, `absolute_cnv`,
            `relative_cnv`, `repeated_seq_expr`, `literal_seq_expr`. This parameter
            determines how to represent HGVS dup/del expressions as VRS objects.
        :param Optionals[Tuple] ival: Interval
        :param bool use_vrs_allele_range: `True` if allele should be computed
            using `to_vrs_allele_ranges` method. `False` if allele should be
            computed using `to_vrs_allele` method.
        :param Optional[int] baseline_copies: Baseline copies number
        :param Optional[RelativeCopyClass] relative_copy_class: The relative copy class
        """
        if errors:
            return

        if grch38:
            t = grch38["ac"]

        if use_vrs_allele_range:
            allele = self.vrs.to_vrs_allele_ranges(
                t, s.coordinate_type, s.alt_type, errors, ival)
        else:
            allele = self.vrs.to_vrs_allele(
                t, grch38["pos"][0], grch38["pos"][1], s.coordinate_type,
                s.alt_type, errors)

        grch38_variation = self.hgvs_dup_del_mode.interpret_variation(
            s.alt_type, allele, errors, hgvs_dup_del_mode,
            baseline_copies=baseline_copies, relative_copy_class=relative_copy_class)

        if grch38_variation:
            self._add_dict_to_mane_data(
                grch38["ac"], s, grch38_variation, mane_data_found, "GRCh38")
