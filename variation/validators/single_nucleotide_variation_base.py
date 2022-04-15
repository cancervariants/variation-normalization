"""The module for Single Nucleotide Variation Validation."""
from typing import List, Dict, Optional
import logging

from ga4gh.vrsatile.pydantic.vrs_models import RelativeCopyClass

from variation.schemas.app_schemas import Endpoint
from variation.schemas.classification_response_schema import Classification, \
    ClassificationType
from variation.schemas.token_response_schema import Token, GeneMatchToken
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum
from .validator import Validator

logger = logging.getLogger("variation")
logger.setLevel(logging.DEBUG)


class SingleNucleotideVariationBase(Validator):
    """The Single Nucleotide Variation Validator Base class."""

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

    def get_gene_tokens(
            self, classification: Classification) -> List[GeneMatchToken]:
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

    async def silent_mutation_valid_invalid_results(
            self, classification_tokens: List, transcripts: List,
            classification: Classification, results: List, gene_tokens: List,
            endpoint_name: Optional[Endpoint], mane_data_found: Dict,
            is_identifier: bool, do_liftover: bool = False) -> None:
        """Add validation result objects to a list of results for
        Silent Mutations.

        :param List classification_tokens: A list of classification Tokens
        :param List transcripts: A list of transcript accessions
        :param Classification classification: A classification for a list of
            tokens
        :param List results: Stores validation result objects
        :param List gene_tokens: List of GeneMatchTokens for a classification
        :param Optional[Endpoint] endpoint_name: Then name of the endpoint being used
        :param Dict mane_data_found: MANE Transcript information found
        :param bool is_identifier: `True` if identifier is given for exact
            location. `False` otherwise.
        :param bool do_liftover: Whether or not to liftover to GRCh38 assembly
        """
        valid_alleles = list()
        for s in classification_tokens:
            for t in transcripts:
                errors = list()

                t = self.get_accession(t, classification)
                allele = None
                if s.coordinate_type == "c":
                    cds_start_end = await self.uta.get_cds_start_end(t)

                    if not cds_start_end:
                        cds_start = None
                        errors.append(
                            f"Unable to find CDS start for accession : {t}"
                        )
                    else:
                        cds_start = cds_start_end[0]
                else:
                    cds_start = None

                if not errors:
                    allele = self.vrs.to_vrs_allele(
                        t, s.position, s.position, s.coordinate_type,
                        s.alt_type, errors, cds_start=cds_start,
                        alt=s.new_nucleotide
                    )

                if not errors:
                    sequence, w = self.seqrepo_access.get_reference_sequence(
                        t, s.position)
                    if sequence is None:
                        errors.append(w)
                    else:
                        if s.ref_nucleotide:
                            if sequence != s.ref_nucleotide:
                                errors.append(
                                    f"Expected {s.coordinate_type} but "
                                    f"found {sequence}")

                if not errors:
                    if endpoint_name == Endpoint.NORMALIZE:
                        mane = await self.mane_transcript.get_mane_transcript(
                            t, s.position, s.coordinate_type, end_pos=s.position,
                            gene=gene_tokens[0].token if gene_tokens else None,
                            try_longest_compatible=True)
                        self.add_mane_data(mane, mane_data_found, s.coordinate_type,
                                           s.alt_type, s)
                    elif endpoint_name == Endpoint.TO_CANONICAL and do_liftover:
                        await self._liftover_genomic_data(
                            gene_tokens, t, s, errors, valid_alleles, results,
                            classification)

                self.add_validation_result(allele, valid_alleles, results,
                                           classification, s, t, gene_tokens, errors)

                if is_identifier:
                    break

        if endpoint_name == Endpoint.NORMALIZE:
            self.add_mane_to_validation_results(mane_data_found, valid_alleles, results,
                                                classification, gene_tokens)

    async def _liftover_genomic_data(
        self, gene_tokens: List, t: str, s: Token, errors: List, valid_alleles: List,
        results: List, classification: Classification
    ) -> None:
        """Add liftover data to validation results
        Currently only works for HGVS expressions.

        :param List gene_tokens: List of GeneMatchTokens for a classification
        :param str t: Accession
        :param Token s: Classification token
        :param List errors: List of errors
        :param List valid_alleles: List of valid alleles
        :param List results: List of results data
        :param Classification classification: A classification for a list of tokens
        """
        if not gene_tokens:
            if not self._is_grch38_assembly(t):
                grch38 = await self.mane_transcript.g_to_grch38(t, s.position,
                                                                s.position)
            else:
                grch38 = dict(ac=t, pos=(s.position, s.position))

            await self.add_genomic_liftover_to_results(
                grch38, errors, s.new_nucleotide, valid_alleles, results,
                classification, s, t, gene_tokens)

    def check_ref_nucleotide(self, actual_ref_nuc: str, expected_ref_nuc: str,
                             position: int, t: str, errors: List) -> None:
        """Assert that ref_nuc matches s.ref_nucleotide."""
        if actual_ref_nuc != expected_ref_nuc:
            errors.append(f"Needed to find {expected_ref_nuc} at position {position} "
                          f"on {t} but found {actual_ref_nuc}")
