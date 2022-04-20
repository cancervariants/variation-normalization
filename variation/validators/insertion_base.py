"""The module for Insertion Validation."""
from typing import List, Dict, Optional
import logging

from ga4gh.vrsatile.pydantic.vrs_models import RelativeCopyClass

from variation.schemas.classification_response_schema import Classification
from variation.schemas.app_schemas import Endpoint
from variation.schemas.token_response_schema import Token
from variation.validators.validator import Validator
from variation.schemas.normalize_response_schema\
    import HGVSDupDelMode as HGVSDupDelModeEnum

logger = logging.getLogger("variation")
logger.setLevel(logging.DEBUG)


class InsertionBase(Validator):
    """The Insertion Validator Base class."""

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
                allele = None

                if s.coordinate_type == "c":
                    cds_start_end = await self.uta.get_cds_start_end(t)
                    if cds_start_end is not None:
                        cds_start = cds_start_end[0]
                    else:
                        cds_start = 0
                        allele = None
                        errors.append(f"Unable to get CDS start for {t}")
                else:
                    cds_start = None

                if not errors:
                    allele = self.vrs.to_vrs_allele(
                        t, s.start_pos_flank, s.end_pos_flank,
                        s.coordinate_type, s.alt_type, errors,
                        cds_start=cds_start, alt=s.inserted_sequence
                    )

                if not errors:
                    self.check_pos_index(t, s, errors)

                if not errors:
                    if endpoint_name == Endpoint.NORMALIZE:
                        mane = await self.mane_transcript.get_mane_transcript(
                            t, s.start_pos_flank, s.coordinate_type,
                            end_pos=s.end_pos_flank,
                            gene=gene_tokens[0].token if gene_tokens else None,
                            try_longest_compatible=True
                        )
                        self.add_mane_data(mane, mane_data_found, s.coordinate_type,
                                           s.alt_type, s, alt=s.inserted_sequence)
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
        """Add liftover data to validation results.
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
                grch38 = await self.mane_transcript.g_to_grch38(t, s.start_pos_flank,
                                                                s.end_pos_flank)
            else:
                grch38 = dict(ac=t, pos=(s.start_pos_flank, s.end_pos_flank))

            await self.add_genomic_liftover_to_results(
                grch38, errors, s.inserted_sequence, valid_alleles, results,
                classification, s, t, gene_tokens)

    def get_hgvs_expr(self, classification: Classification, t: str, s: Token,
                      is_hgvs: bool) -> str:
        """Return a HGVS expression.

        :param Classification classification: A classification for a list of
            tokens
        :param str t: Transcript retrieved from transcript mapping
        :param Token s: The classification token
        :param bool is_hgvs: Whether or not classification is HGVS token
        :return: hgvs expression
        """
        if not is_hgvs:
            prefix = f"{t}:{s.coordinate_type.lower()}."
            position = f"{s.start_pos_flank}_{s.end_pos_flank}"
            if s.inserted_sequence2 is not None:
                inserted_sequence = f"{s.inserted_sequence}_{s.inserted_sequence2}"
            else:
                inserted_sequence = f"{s.inserted_sequence}"

            hgvs_expr = f"{prefix}{position}ins{inserted_sequence}"
        else:
            hgvs_token = [t for t in classification.all_tokens if
                          isinstance(t, Token) and t.token_type
                          in ["HGVS", "ReferenceSequence"]][0]
            hgvs_expr = hgvs_token.input_string
        return hgvs_expr

    def check_pos_index(self, t: str, s: Token, errors: List) -> None:
        """Check that position exists on transcript.

        :param str t: Transcript accession
        :param Token s: Classification token
        :param List errors: List of errors
        """
        sequence, w = self.seqrepo_access.get_reference_sequence(
            t, int(s.start_pos_flank), int(s.end_pos_flank))

        if sequence is None:
            errors.append(w)
