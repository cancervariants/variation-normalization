"""The module for Protein Insertion Validation."""
import logging
from typing import List, Optional, Dict

from ga4gh.vrsatile.pydantic.vrs_models import RelativeCopyClass
from ga4gh.vrs.extras.translator import Translator
from gene.query import QueryHandler as GeneQueryHandler
from uta_tools.data_sources import SeqRepoAccess, TranscriptMappings, UTADatabase, \
    MANETranscript

from variation.schemas.classification_response_schema import \
    ClassificationType, Classification
from variation.schemas.token_response_schema import Token
from variation.schemas.app_schemas import Endpoint
from variation.validators.validator import Validator
from variation.schemas.token_response_schema import GeneMatchToken
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import AminoAcidCache
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum
from variation.vrs_representation import VRSRepresentation
from .protein_base import ProteinBase


logger = logging.getLogger("variation")
logger.setLevel(logging.DEBUG)


class ProteinInsertion(Validator):
    """The Protein Insertion Validator class."""

    def __init__(
        self, seq_repo_access: SeqRepoAccess, transcript_mappings: TranscriptMappings,
        gene_symbol: GeneSymbol, mane_transcript: MANETranscript, uta: UTADatabase,
        tlr: Translator, gene_normalizer: GeneQueryHandler,
        vrs: VRSRepresentation, amino_acid_cache: AminoAcidCache
    ) -> None:
        """Initialize the validator.

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
        :param AminoAcidCache amino_acid_cache: Amino Acid codes and conversions
        """
        super().__init__(
            seq_repo_access, transcript_mappings, gene_symbol, mane_transcript,
            uta, tlr, gene_normalizer, vrs
        )
        self._amino_acid_cache = amino_acid_cache
        self.protein_base = ProteinBase(seq_repo_access, amino_acid_cache)
        self.mane_transcript = mane_transcript

    async def get_transcripts(self, gene_tokens: List, classification: Classification,
                              errors: List) -> Optional[List[str]]:
        """Get transcript accessions for a given classification.

        :param List gene_tokens: A list of gene tokens
        :param Classification classification: A classification for a list of
            tokens
        :param List errors: List of errors
        :return: List of transcript accessions
        """
        return self.get_protein_transcripts(gene_tokens, errors)

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
        valid_alleles = list()
        for s in classification_tokens:
            for t in transcripts:
                errors = list()
                t = self.get_accession(t, classification)
                allele = self.vrs.to_vrs_allele(
                    t, s.start_pos_flank, s.end_pos_flank,
                    s.coordinate_type, s.alt_type, errors,
                    alt=s.inserted_sequence
                )

                if not errors:
                    # Check ref start/end protein matches expected
                    self.protein_base.check_ref_aa(
                        t, s.start_aa_flank, s.start_pos_flank, errors
                    )

                    if not errors:
                        self.protein_base.check_ref_aa(
                            t, s.end_aa_flank, s.end_pos_flank, errors
                        )

                if not errors and endpoint_name == Endpoint.NORMALIZE:
                    mane = await self.mane_transcript.get_mane_transcript(
                        t, s.start_pos_flank, s.coordinate_type,
                        end_pos=s.end_pos_flank, try_longest_compatible=True)
                    self.add_mane_data(mane, mane_data_found, s.coordinate_type,
                                       s.alt_type, s, alt=s.inserted_sequence)

                self.add_validation_result(allele, valid_alleles, results,
                                           classification, s, t, gene_tokens, errors)

                if is_identifier:
                    break

        if endpoint_name == Endpoint.NORMALIZE:
            self.add_mane_to_validation_results(mane_data_found, valid_alleles, results,
                                                classification, gene_tokens)

    def get_gene_tokens(self, classification: Classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_protein_gene_symbol_tokens(classification)

    def variation_name(self) -> str:
        """Return the variation name."""
        return "protein insertion"

    def is_token_instance(self, t: Token) -> bool:
        """Check that token is Protein Insertion."""
        return t.token_type == "ProteinInsertion"

    def validates_classification_type(
        self, classification_type: ClassificationType
    ) -> bool:
        """Return whether or not the classification type is
        Protein Insertion.
        """
        return classification_type == ClassificationType.PROTEIN_INSERTION
