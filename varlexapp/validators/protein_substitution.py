"""The module for Protein Substitution Validation."""
from typing import List
from requests.exceptions import HTTPError
from .validator import Validator
from varlexapp.schemas.classification_response_schema import \
    ClassificationType, Classification
from varlexapp.schemas.token_response_schema import ProteinSubstitutionToken
from varlexapp.schemas.validation_response_schema import ValidationResult, \
    LookupType
from varlexapp.schemas.token_response_schema import Token
from varlexapp.tokenizers import GeneSymbol
from varlexapp.tokenizers.caches import GeneSymbolCache
from varlexapp.tokenizers.caches import AminoAcidCache
from varlexapp.data_sources import SeqRepoAccess, TranscriptMappings
from ga4gh.vrs import models
from ga4gh.core import ga4gh_identify
from ga4gh.vrs.dataproxy import SeqRepoRESTDataProxy
from ga4gh.vrs.extras.translator import Translator


class ProteinSubstitution(Validator):
    """The Protein Substitution Validator class."""

    def __init__(self, seq_repo_client: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 seqrepo_rest_service_url="http://localhost:5000/seqrepo") \
            -> None:
        """Initialize the Protein Substitution validator."""
        self.transcript_mappings = transcript_mappings
        self.seq_repo_client = seq_repo_client
        self._gene_matcher = GeneSymbol(GeneSymbolCache())
        self._amino_acid_cache = AminoAcidCache()
        self.seqrepo_rest_service = seqrepo_rest_service_url
        self.dp = SeqRepoRESTDataProxy(base_url=seqrepo_rest_service_url)
        self.tlr = Translator(data_proxy=self.dp)

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
            return [ValidationResult(
                classification=classification,
                is_valid=False,
                confidence_score=0,
                allele=None,
                human_description='',
                concise_description='',
                errors=errors
            )]

        transcripts = self.transcript_mappings.protein_transcripts(
            gene_tokens[0].token, LookupType.GENE_SYMBOL)

        if not transcripts:
            errors.append(f'No transcripts found for gene symbol '
                          f'{gene_tokens[0].token}')
            return [ValidationResult(
                classification=classification,
                is_valid=False,
                confidence_score=0,
                allele=None,
                human_description='',
                concise_description='',
                errors=errors
            )]

        for s in psub_tokens:
            for t in transcripts:
                valid = True
                errors = list()
                ref_protein = \
                    self.seq_repo_client.protein_at_position(t, s.position)

                # TODO: Should sequence_id be GS or SQ?

                if 'HGVS' in classification.matching_tokens:
                    hgvs_token = \
                        [t for t in classification.all_tokens if
                         isinstance(t, Token) and t.token_type == 'HGVS'][0]
                    hgvs_expr = hgvs_token.input_string

                    try:
                        seq_location = self.tlr.translate_from(hgvs_expr,
                                                               'hgvs')
                    except HTTPError:
                        errors.append(f"{hgvs_expr} is an invalid "
                                      f"HGVS expression.")
                        valid = False
                        allele = None
                    else:
                        allele = seq_location.as_dict()
                        allele['location']['sequence_id'] =\
                            allele['location']['sequence_id'].replace(
                                "ga4gh:GS", "ga4gh:SQ")
                else:
                    sequence_id = \
                        self.dp.translate_sequence_identifier(
                            t, 'ga4gh')[0].replace("ga4gh:GS", "ga4gh:SQ")

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
                    allele = allele.as_dict()

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
                    results.append(ValidationResult(
                        classification=classification,
                        is_valid=True,
                        confidence_score=1,
                        allele=allele,
                        human_description=self.concise_description(t, s),
                        concise_description=self.concise_description(t, s),
                        errors=[]
                    ))
                else:
                    results.append(ValidationResult(
                        classification=classification,
                        is_valid=False,
                        confidence_score=1,
                        allele=allele,
                        human_description=self.concise_description(t, s),
                        concise_description=self.concise_description(t, s),
                        errors=errors
                    ))
                errors = list()
        return results

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
