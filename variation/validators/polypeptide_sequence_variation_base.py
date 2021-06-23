"""The module for Polypeptide Sequence Variation Validation."""
from typing import List, Optional
from abc import abstractmethod
from .validator import Validator
from variation.schemas.token_response_schema import GeneMatchToken
from variation.schemas.token_response_schema import Token
from variation.tokenizers import GeneSymbol
from variation.tokenizers.caches import AminoAcidCache
from variation.data_sources import SeqRepoAccess, TranscriptMappings
from variation.mane_transcript import MANETranscript
from .amino_acid_base import AminoAcidBase
import logging

logger = logging.getLogger('variation')
logger.setLevel(logging.DEBUG)


class PolypeptideSequenceVariationBase(Validator):
    """The Polypeptide Sequence Variation Validator Base class."""

    def __init__(self, seq_repo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol,
                 mane_transcript: MANETranscript,
                 amino_acid_cache: AminoAcidCache) \
            -> None:
        """Initialize the validator.

        :param SeqRepoAccess seq_repo_access: Access to SeqRepo data
        :param TranscriptMappings transcript_mappings: Access to transcript
            mappings
        """
        super().__init__(
            seq_repo_access, transcript_mappings, gene_symbol, mane_transcript
        )
        self._amino_acid_cache = amino_acid_cache
        self.amino_acid_base = AminoAcidBase(seq_repo_access, amino_acid_cache)
        self.mane_transcript = mane_transcript

    def get_transcripts(self, gene_tokens, classification, errors)\
            -> Optional[List[str]]:
        """Get transcript accessions for a given classification.

        :param list gene_tokens: A list of gene tokens
        :param Classification classification: A classification for a list of
            tokens
        :param list errors: List of errors
        :return: List of transcript accessions
        """
        return self.get_protein_transcripts(gene_tokens, errors)

    def get_hgvs_expr(self, classification, t, s, is_hgvs) -> str:
        """Return HGVS expression for a classification.

        :param Classification classification: A classification for a list of
            tokens
        :param str t: Transcript retrieved from transcript mapping
        :param Token s: The classification token
        :param bool is_hgvs: Whether or not classification is HGVS token
        :return: hgvs expression
        """
        if not is_hgvs:
            hgvs_expr = f"{t}:p.{s.ref_protein}{s.position}" \
                        f"{s.alt_protein}"
        else:
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
                    self._amino_acid_cache.amino_acid_code_conversion[alt_amino]  # noqa: E501
                hgvs_expr = hgvs_expr.replace('=', three_letter)
        return hgvs_expr

    def get_valid_invalid_results(self, classification_tokens, transcripts,
                                  classification, results, gene_tokens,
                                  normalize_endpoint) -> None:
        """Add validation result objects to a list of results.

        :param list classification_tokens: A list of Tokens
        :param list transcripts: A list of transcript strings
        :param Classification classification: A classification for a list of
            tokens
        :param list results: A list to store validation result objects
        :param list gene_tokens: List of GeneMatchTokens
        """
        valid_alleles = list()
        if 'HGVS' in classification.matching_tokens:
            is_hgvs = True
        else:
            is_hgvs = False

        mane_data = {
            'mane_select': list(),
            'mane_plus_clinical': list(),
            'longest_compatible_remaining': list()
        }

        for s in classification_tokens:
            for t in transcripts:
                errors = list()
                allele, t, hgvs_expr, is_ensembl = \
                    self.get_allele_with_context(classification, t, s, errors)

                mane = self.mane_transcript.get_mane_transcript(
                    t, s.position, s.position, 'p',
                    ref=self._amino_acid_cache.convert_three_to_one(
                        s.ref_protein
                    ),
                    normalize_endpoint=normalize_endpoint
                )
                if mane:
                    hgvs_expr = f"{mane['refseq']}:p.{s.ref_protein}" \
                                f"{mane['pos'][0]}{s.alt_protein}"
                    key = '_'.join(mane['status'].lower().split())
                    mane_data[key].append(
                        {
                            'hgvs_expr': hgvs_expr,
                            'classification_token': s,
                            'transcript': mane['refseq']
                        }
                    )

                if not allele:
                    errors.append("Unable to find allele.")
                else:
                    if len(allele['state']['sequence']) == 3:
                        allele['state']['sequence'] = \
                            self._amino_acid_cache.convert_three_to_one(
                                allele['state']['sequence'])

                    self.amino_acid_base.check_ref_aa(
                        t, s.ref_protein, s.position, errors
                    )

                self.add_validation_result(allele, valid_alleles, results,
                                           classification, s, t, gene_tokens,
                                           errors)
                if is_hgvs:
                    break

        for key in ['mane_select', 'mane_plus_clinical',
                    'longest_compatible_remaining']:
            for data in mane_data[key][::-1]:
                errors = list()
                found_mane = self.add_validation_result(
                    self.get_allele_from_hgvs(data['hgvs_expr'], errors),
                    valid_alleles, results, classification,
                    data['classification_token'], data['transcript'],
                    gene_tokens, errors, mane_transcript=data['hgvs_expr']
                )
                if found_mane:
                    return

    def get_gene_tokens(self, classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_protein_gene_symbol_tokens(classification)

    def concise_description(self, transcript, token) -> str:
        """Return a HGVS description of the identified variation.

        :param str transcript: Transcript accession
        :param Token token: Classification token
        :return: HGVS expression
        """
        return f'{transcript} {token.ref_protein}' \
               f'{token.position}{token.alt_protein}'

    @abstractmethod
    def human_description(self, transcript, token) -> str:
        """Return a human description of the identified variation."""
        raise NotImplementedError
