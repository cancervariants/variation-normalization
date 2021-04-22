"""The module for Genomic Substitution Validation."""
from .single_nucleotide_variant_base import SingleNucleotideVariantBase
from variant.schemas.classification_response_schema import \
    ClassificationType
from variant.schemas.token_response_schema import GenomicSubstitutionToken
from typing import List
from variant.schemas.classification_response_schema import Classification
from variant.schemas.validation_response_schema import ValidationResult
import logging
from variant.schemas.token_response_schema import Token
from gene.query import QueryHandler as GeneQueryHandler

logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)

# TODO: Find gene from NC accession (in event of no mane transcripts)


class GenomicSubstitution(SingleNucleotideVariantBase):
    """The Genomic Substitution Validator class."""

    def validate(self, classification: Classification) \
            -> List[ValidationResult]:
        """Validate a given classification.

        :param Classification classification: A classification for a list of
            tokens
        :return: A list of validation results
        """
        results = list()
        errors = list()

        classification_tokens = self.get_classification_tokens(classification)
        gene_tokens = self.get_gene_tokens(classification)

        if gene_tokens and len(gene_tokens) > 1:
            errors.append('More than one gene symbol found for a single'
                          f' {self.variant_name()}')

        if len(classification.non_matching_tokens) > 0:
            errors.append(f"Non matching tokens found for "
                          f"{self.variant_name()}.")

        nc_accessions = self.get_nc_accessions(classification)
        if not nc_accessions:
            errors.append('Could not find NC_ accession for '
                          f'{self.variant_name()}')

        if len(errors) > 0:
            return [self.get_validation_result(
                classification, False, 0, None,
                '', '', errors, gene_tokens)]

        self.get_valid_invalid_results(classification_tokens, nc_accessions,
                                       classification, results, gene_tokens)
        return results

    def get_hgvs_expr(self, classification, t):
        """Get HGVS expression."""
        hgvs_token = [t for t in classification.all_tokens if
                      isinstance(t, Token) and t.token_type == 'HGVS'][0]
        input_hgvs_expr = hgvs_token.input_string.split(':')[0]
        if input_hgvs_expr != t:
            hgvs_token = f"{t}:{hgvs_token.input_string.split(':')[1]}"
        else:
            hgvs_token = hgvs_token.input_string
        return hgvs_token

    def get_valid_invalid_results(self, classification_tokens, transcripts,
                                  classification, results,
                                  gene_tokens) -> None:
        """Add validation result objects to a list of results.

        :param list classification_tokens: A list of Tokens
        :param list transcripts: A list of transcript strings
        :param Classification classification: A classification for a list of
            tokens
        :param list results: A list to store validation result objects
        :param list gene_tokens: List of GeneMatchTokens
        """
        valid_alleles = list()
        mane_transcripts_dict = dict()
        for s in classification_tokens:
            for t in transcripts:
                errors = list()
                ref_nuc = \
                    self.seq_repo_access.sequence_at_position(t, s.position)

                if 'HGVS' in classification.matching_tokens:
                    hgvs_expr = self.get_hgvs_expr(classification, t)
                    is_ensembl_transcript = False
                    allele = self.get_allele_from_hgvs(hgvs_expr, errors)
                    if allele:
                        mane_transcripts_dict[hgvs_expr] = {
                            'classification_token': s,
                            'transcript_token': t,
                            'is_ensembl_transcript': is_ensembl_transcript
                        }
                else:
                    allele = self.get_allele_from_transcript(s, t, errors)
                    if allele:
                        self._add_hgvs_to_mane_transcripts_dict(
                            classification, mane_transcripts_dict, s, t,
                            gene_tokens
                        )
                self.check_ref_nucleotide(ref_nuc, s, t, errors)
                self.add_validation_result(
                    allele, valid_alleles, results,
                    classification, s, t, gene_tokens, errors
                )

        # Now add Mane transcripts to results
        self.add_mane_transcript(classification, results, gene_tokens,
                                 mane_transcripts_dict)

    def get_gene_tokens(self, classification):
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_gene_symbol_tokens(classification)

    def get_nc_accessions(self, classification):
        """Get NC accession for a given classification."""
        hgvs = [t.token for t in classification.all_tokens if
                t.token_type in ['HGVS', 'ReferenceSequence']]
        nc_accessions = []
        if hgvs:
            nc_accession = hgvs[0].split(':')[0]
            nc_accessions = \
                self.get_nc_accessions_from_nc_accession(nc_accession)
        else:
            gene_tokens = [t for t in classification.all_tokens
                           if t.token_type == 'GeneSymbol']
            if gene_tokens and len(gene_tokens) == 1:
                gene_query_handler = GeneQueryHandler()
                resp = gene_query_handler.search_sources(gene_tokens[0].token,
                                                         incl='hgnc')
                if resp['source_matches'][0]['records']:
                    record = resp['source_matches'][0]['records'][0]
                    loc = record.locations[0] if record.locations else None
                    # TODO: what about multiple chr locations?
                    if loc and loc.chr:
                        for identifier in ['GRCh38', 'GRCh37']:
                            nc_accession =  \
                                self.get_nc_accession(f"{identifier}:"
                                                      f"{loc.chr}")
                            if nc_accession:
                                nc_accessions.append(nc_accession)
        return list(set(nc_accessions))

    def get_nc_accessions_from_nc_accession(self, nc_accession):
        """Given NC accession, find other version from other assembly."""
        nc_accessions = [nc_accession]
        try:
            assembly = None
            for a in self.dp.get_metadata(nc_accession)['aliases']:
                if a.startswith('GRCh3'):
                    assembly = a
                    break
        except KeyError:
            pass
        else:
            if assembly:
                if assembly.startswith('GRCh38'):
                    nc_accession = \
                        self.get_nc_accession(f"GRCh37:"
                                              f"{assembly.split(':')[1]}")
                    if nc_accession:
                        nc_accessions.append(nc_accession)
                elif assembly.startswith('GRCh37'):
                    nc_accession = self.get_nc_accession(
                        f"GRCh38:{assembly.split(':')[1]}")
                    if nc_accession:
                        nc_accessions.append(nc_accession)
        return nc_accessions

    def get_nc_accession(self, identifier):
        """Given an identifier (assembly+chr), return nc accession."""
        nc_accession = None
        try:
            metadata = \
                self.dp.get_metadata(identifier)
        except KeyError:
            logger.warning('Data Proxy unable to get metadata'
                           f'for GRCh38:{identifier}')
        else:
            aliases = [a for a in metadata['aliases'] if
                       a.startswith('refseq:NC_')]
            if aliases:
                nc_accession = aliases[0].split(':')[-1]
        return nc_accession

    def variant_name(self):
        """Return the variant name."""
        return 'genomic substitution'

    def is_token_instance(self, t):
        """Check that token is Coding DNA Substitution."""
        return t.token_type == 'GenomicSubstitution'

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return whether or not the classification type is amino acid
        substitution.
        """
        return classification_type == ClassificationType.GENOMIC_SUBSTITUTION

    def human_description(self, transcript,
                          psub_token: GenomicSubstitutionToken) -> str:
        """Return a human description of the identified variant."""
        return f'A genomic DNA substitution from {psub_token.ref_nucleotide}' \
               f' to {psub_token.new_nucleotide} at position ' \
               f'{psub_token.position} on transcript {transcript}'