"""The module for Protein Substitution Validation."""
from typing import List
from requests.exceptions import HTTPError
from .validator import Validator
from variant.schemas.classification_response_schema import \
    ClassificationType, Classification
from variant.schemas.token_response_schema import ProteinSubstitutionToken
from variant.schemas.validation_response_schema import ValidationResult, \
    LookupType
from variant.schemas.token_response_schema import Token
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import GeneSymbolCache, AminoAcidCache
from variant.data_sources import SeqRepoAccess, TranscriptMappings
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify
from ga4gh.vrs.dataproxy import SeqRepoDataProxy
from ga4gh.vrs.extras.translator import Translator
from variant import SEQREPO_DATA_PATH
from biocommons.seqrepo import SeqRepo
import hgvs.parser


class ProteinSubstitution(Validator):
    """The Protein Substitution Validator class."""

    def __init__(self, seq_repo_client: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings) \
            -> None:
        """Initialize the Protein Substitution validator."""
        self.transcript_mappings = transcript_mappings
        self.seq_repo_client = seq_repo_client
        self._gene_matcher = GeneSymbol(GeneSymbolCache())
        self._amino_acid_cache = AminoAcidCache()
        self.dp = SeqRepoDataProxy(SeqRepo(SEQREPO_DATA_PATH))
        self.tlr = Translator(data_proxy=self.dp)
        self.hgvs_parser = hgvs.parser.Parser()

    def validate(self, classification: Classification) \
            -> List[ValidationResult]:
        """Validate a given protein substitution classification."""
        results = list()
        errors = list()

        psub_tokens = [t for t in classification.all_tokens if isinstance(t, ProteinSubstitutionToken)]  # noqa: E501

        gene_tokens = self.get_gene_tokens(classification)

        if len(gene_tokens) == 0:
            errors.append('No gene tokens for a protein substitution.')

        if len(gene_tokens) > 1:
            errors.append('More than one gene symbol found for a single'
                          ' protein substitution')

        if len(errors) > 0:
            return [self.get_validation_result(
                        classification, False, 0, None,
                        '', '', errors)]

        transcripts = self.transcript_mappings.protein_transcripts(
            gene_tokens[0].token, LookupType.GENE_SYMBOL)

        if not transcripts:
            errors.append(f'No transcripts found for gene symbol '
                          f'{gene_tokens[0].token}')
            return [self.get_validation_result(
                        classification, False, 0, None,
                        '', '', errors)]

        self.get_valid_invalid_results(psub_tokens, transcripts,
                                       classification, results)
        return results

    def get_vrs_allele(self, sequence_id, s):
        """Return VRS Allele object."""
        seq_location = models.SequenceLocation(
            sequence_id=sequence_id,
            interval=models.SimpleInterval(
                start=s.position - 1,
                end=s.position
            )
        )

        state = models.SequenceState(sequence=s.alt_protein)
        allele = models.Allele(location=seq_location,
                               state=state)
        allele['_id'] = ga4gh_identify(allele)
        return allele.as_dict()

    def get_hgvs_expr(self, classification):
        """Return HGVS expression."""
        # Replace `=` in silent mutation with 3 letter amino acid code
        hgvs_token = \
            [t for t in classification.all_tokens if
             isinstance(t, Token) and t.token_type == 'HGVS'][0]
        hgvs_expr = hgvs_token.input_string

        if '=' in hgvs_expr:
            hgvs_parsed = \
                self.hgvs_parser.parse_hgvs_variant(hgvs_expr)
            alt_amino = hgvs_parsed.posedit.pos.end.aa
            three_letter = \
                self._amino_acid_cache._amino_acid_code_conversion[alt_amino]
            hgvs_expr = hgvs_expr.replace('=', three_letter)
        return hgvs_expr

    def get_valid_invalid_results(self, psub_tokens, transcripts,
                                  classification, results):
        """Add validation result objects to a list of results."""
        for s in psub_tokens:
            for t in transcripts:
                valid = True
                errors = list()
                ref_protein = \
                    self.seq_repo_client.protein_at_position(t, s.position)

                if 'HGVS' in classification.matching_tokens:
                    hgvs_expr = self.get_hgvs_expr(classification)

                    try:
                        # TODO: We should get this to work w/o doing above
                        #       replacement
                        seq_location = self.tlr.translate_from(hgvs_expr,
                                                               'hgvs')
                    except HTTPError:
                        errors.append(f"{hgvs_expr} is an invalid "
                                      f"HGVS expression.")
                        valid = False
                        allele = None
                    else:
                        allele = seq_location.as_dict()
                else:
                    sequence_id = \
                        self.dp.translate_sequence_identifier(t, 'ga4gh')[0]
                    allele = self.get_vrs_allele(sequence_id, s)

                if not errors:
                    if ref_protein and len(ref_protein) == 1 \
                            and len(s.ref_protein) == 3:
                        ref_protein = self._amino_acid_cache._amino_acid_code_conversion[ref_protein]  # noqa: E501
                    if ref_protein != s.ref_protein:
                        errors.append(f'Needed to find {s.ref_protein} at'
                                      f' position {s.position} on {t}'
                                      f' but found {ref_protein}')
                        valid = False

                if valid:
                    results.append(self.get_validation_result(
                        classification, True, 1, allele,
                        self.human_description(t, s),
                        self.concise_description(t, s), []))
                else:
                    results.append(self.get_validation_result(
                        classification, False, 1, allele,
                        self.human_description(t, s),
                        self.concise_description(t, s), errors))
                errors = list()

    def get_validation_result(self, classification, is_valid, confidence_score,
                              allele, human_description, concise_description,
                              errors):
        """Return a validation result object."""
        return ValidationResult(
            classification=classification,
            is_valid=is_valid,
            confidence_score=confidence_score,
            allele=allele,
            human_description=human_description,
            concise_description=concise_description,
            errors=errors
        )

    def get_gene_tokens(self, classification):
        """Return gene tokens for a classification."""
        gene_tokens = [t for t in classification.all_tokens
                       if t.token_type == 'GeneSymbol']
        if not gene_tokens:
            # Convert refseq to gene symbol
            refseq = \
                [t.token for t in classification.all_tokens if
                 t.token_type in ['HGVS', 'ReferenceSequence']][0]

            if ':' in refseq:
                refseq = refseq.split(':')[0]

            res = self.seq_repo_client.aliases(refseq)
            aliases = [a.split('ensembl:')[1] for a
                       in res if a.startswith('ensembl')]

            gene_symbols = list()
            for alias in aliases:
                gene_symbol = \
                    self.transcript_mappings.refseq_gene_symbol(alias)

                if gene_symbol and gene_symbol not in gene_symbols:
                    gene_symbols.append(gene_symbol)
                    gene_tokens.append(self._gene_matcher.match(gene_symbol))
        return gene_tokens

    def validates_classification_type(
            self,
            classification_type: ClassificationType) -> bool:
        """Return the protein substitution classification type."""
        return classification_type == ClassificationType.PROTEIN_SUBSTITUTION

    def concise_description(self, transcript,
                            psub_token: ProteinSubstitutionToken) -> str:
        """Return a precise description."""
        return f'{transcript} {psub_token.ref_protein}' \
               f'{psub_token.position}{psub_token.alt_protein}'

    def human_description(self, transcript,
                          psub_token: ProteinSubstitutionToken) -> str:
        """Return a human description."""
        return f'A protein substitution from {psub_token.ref_protein}' \
               f' to {psub_token.alt_protein} at position ' \
               f'{psub_token.position} on transcript {transcript}'
