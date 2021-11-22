"""The module for Single Nucleotide Variation Validation."""
from typing import List, Dict
from .validator import Validator
import logging
from variation.schemas.classification_response_schema import Classification

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class SingleNucleotideVariationBase(Validator):
    """The Single Nucleotide Variation Validator Base class."""

    def silent_mutation_valid_invalid_results(
            self, classification_tokens: List, transcripts: List,
            classification: Classification, results: List, gene_tokens: List,
            normalize_endpoint: bool, mane_data_found: Dict,
            is_identifier: bool) -> None:
        """Add validation result objects to a list of results for
        Silent Mutations.

        :param List classification_tokens: A list of classification Tokens
        :param List transcripts: A list of transcript accessions
        :param Classification classification: A classification for a list of
            tokens
        :param List results: Stores validation result objects
        :param List gene_tokens: List of GeneMatchTokens for a classification
        :param bool normalize_endpoint: `True` if normalize endpoint is being
            used. `False` otherwise.
        :param Dict mane_data_found: MANE Transcript information found
        :param bool is_identifier: `True` if identifier is given for exact
            location. `False` otherwise.
        """
        valid_alleles = list()
        for s in classification_tokens:
            for t in transcripts:
                errors = list()

                t = self.get_accession(t, classification)
                allele = None
                if s.reference_sequence == 'c':
                    cds_start_end = self.uta.get_cds_start_end(t)

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
                    allele = self.to_vrs_allele(
                        t, s.position, s.position, s.reference_sequence,
                        s.alt_type, errors, cds_start=cds_start,
                        alt=s.new_nucleotide
                    )

                if not errors:
                    sequence = \
                        self.seqrepo_access.get_sequence(t, s.position)
                    if sequence is None:
                        errors.append('Sequence index error')

                if not errors:
                    mane = self.mane_transcript.get_mane_transcript(
                        t, s.position, s.position, s.reference_sequence,
                        gene=gene_tokens[0].token if gene_tokens else None,
                        normalize_endpoint=normalize_endpoint
                    )

                    self.add_mane_data(mane, mane_data_found,
                                       s.reference_sequence, s.alt_type, s)

                self.add_validation_result(
                    allele, valid_alleles, results,
                    classification, s, t, gene_tokens, errors
                )

                if is_identifier:
                    break
        self.add_mane_to_validation_results(
            mane_data_found, valid_alleles, results,
            classification, gene_tokens
        )

    def check_ref_nucleotide(self, actual_ref_nuc, expected_ref_nuc,
                             position, t, errors):
        """Assert that ref_nuc matches s.ref_nucleotide."""
        if actual_ref_nuc != expected_ref_nuc:
            errors.append(f'Needed to find {expected_ref_nuc} at'
                          f' position {position} on {t}'
                          f' but found {actual_ref_nuc}')

    def concise_description(self, transcript, token) -> str:
        """Return a HGVS description of the identified variation.

        :param str transcript: Transcript accession
        :param Token token: Classification token
        :return: HGVS expression
        """
        prefix = f'{transcript}:{token.reference_sequence}.{token.position}'
        if token.new_nucleotide == '=':
            change = "="
        else:
            change = f"{token.ref_nucleotide}>{token.new_nucleotide}"
        return prefix + change
