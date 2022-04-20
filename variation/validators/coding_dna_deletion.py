"""The module for Coding DNA Deletion Validation."""
import logging
from typing import List, Optional, Dict

from ga4gh.vrsatile.pydantic.vrs_models import RelativeCopyClass

from variation.schemas.app_schemas import Endpoint
from variation.validators.duplication_deletion_base import\
    DuplicationDeletionBase
from variation.schemas.classification_response_schema import \
    ClassificationType, Classification
from variation.schemas.token_response_schema import GeneMatchToken, Token
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum


logger = logging.getLogger("variation")
logger.setLevel(logging.DEBUG)


class CodingDNADeletion(DuplicationDeletionBase):
    """The Coding DNA Deletion Validator class."""

    async def get_transcripts(self, gene_tokens: List, classification: Classification,
                              errors: List) -> Optional[List[str]]:
        """Get transcript accessions for a given classification.

        :param List gene_tokens: A list of gene tokens
        :param Classification classification: A classification for a list of
            tokens
        :param List errors: List of errors
        :return: List of transcript accessions
        """
        return self.get_coding_dna_transcripts(gene_tokens, errors)

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

                cds_start_end = await self.uta.get_cds_start_end(t)
                if cds_start_end is not None:
                    cds_start = cds_start_end[0]

                    allele = self.vrs.to_vrs_allele(
                        t, s.start_pos_del, s.end_pos_del,
                        s.coordinate_type, s.alt_type,
                        errors, cds_start=cds_start
                    )
                else:
                    cds_start = 0
                    allele = None
                    errors.append(f"Unable to find CDS start for {t}")

                if not errors:
                    self.check_reference_sequence(t, s, errors,
                                                  cds_start=cds_start)

                if not errors and endpoint_name == Endpoint.NORMALIZE:
                    mane = await self.mane_transcript.get_mane_transcript(
                        t, s.start_pos_del, s.coordinate_type,
                        end_pos=s.end_pos_del if s.end_pos_del else None,
                        ref=s.deleted_sequence, try_longest_compatible=True,
                        residue_mode="residue")

                    self.add_mane_data(mane, mane_data_found, s.coordinate_type,
                                       s.alt_type, s)
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
        return self.get_coding_dna_gene_symbol_tokens(classification)

    def variation_name(self) -> str:
        """Return the variation name."""
        return "coding dna deletion"

    def is_token_instance(self, t: Token) -> bool:
        """Check that token is Coding DNA Deletion."""
        return t.token_type == "CodingDNADeletion"

    def validates_classification_type(
            self, classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is
        Coding DNA Deletion.
        """
        return classification_type == ClassificationType.CODING_DNA_DELETION
