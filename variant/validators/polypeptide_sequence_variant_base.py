"""The module for Polypeptide Sequence Variant Validation."""
from typing import List, Tuple, Optional
from abc import abstractmethod
from .validator import Validator
from variant.schemas.token_response_schema import GeneMatchToken
from variant.schemas.token_response_schema import Token
from variant.tokenizers import GeneSymbol
from variant.tokenizers.caches import AminoAcidCache
from variant.data_sources import SeqRepoAccess, TranscriptMappings
from .amino_acid_base import AminoAcidBase
import logging

logger = logging.getLogger('variant')
logger.setLevel(logging.DEBUG)


class PolypeptideSequenceVariantBase(Validator):
    """The Polypeptide Sequence Variant Validator Base class."""

    def __init__(self, seq_repo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 gene_symbol: GeneSymbol,
                 amino_acid_cache: AminoAcidCache) \
            -> None:
        """Initialize the validator.

        :param SeqRepoAccess seq_repo_access: Access to SeqRepo data
        :param TranscriptMappings transcript_mappings: Access to transcript
            mappings
        """
        super().__init__(seq_repo_access, transcript_mappings, gene_symbol)
        self._amino_acid_cache = amino_acid_cache
        self.amino_acid_base = AminoAcidBase(seq_repo_access, amino_acid_cache)

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

    def get_hgvs_expr(self, classification, t, s, is_hgvs) -> Tuple[str, None]:
        """Return HGVS expression for a classification.

        :param Classification classification: A classification for a list of
            tokens
        :return: The classification's HGVS expression as a string
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
        return hgvs_expr, None

    def get_valid_invalid_results(self, classification_tokens, transcripts,
                                  classification, results, gene_tokens) \
            -> None:
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

                allele, t, hgvs_expr, is_ensembl_transcript = \
                    self.get_allele_with_context(classification, t, s, errors)

                mane_transcripts_dict[hgvs_expr] = {
                    'classification_token': s,
                    'transcript_token': t
                }

                if allele and len(allele['state']['sequence']) == 3:
                    allele['state']['sequence'] = \
                        self._amino_acid_cache.convert_three_to_one(
                            allele['state']['sequence'])

                if not errors:
                    self.amino_acid_base.check_ref_aa(
                        t, s.ref_protein, s.position, errors
                    )

                if not errors and allele not in valid_alleles:
                    results.append(self.get_validation_result(
                        classification, True, 1, allele,
                        self.human_description(t, s),
                        self.concise_description(t, s), [], gene_tokens))

                    valid_alleles.append(allele)
                else:
                    results.append(self.get_validation_result(
                        classification, False, 1, allele,
                        self.human_description(t, s),
                        self.concise_description(t, s), errors, gene_tokens))

        # Now add Mane transcripts to results
        self.add_mane_transcript(classification, results, gene_tokens,
                                 mane_transcripts_dict)

    def get_gene_tokens(self, classification) -> List[GeneMatchToken]:
        """Return gene tokens for a classification.

        :param Classification classification: The classification for tokens
        :return: A list of Gene Match Tokens in the classification
        """
        return self.get_protein_gene_symbol_tokens(classification)

    def concise_description(self, transcript, token) -> str:
        """Return a description of the identified variant."""
        return f'{transcript} {token.ref_protein}' \
               f'{token.position}{token.alt_protein}'

    @abstractmethod
    def human_description(self, transcript, token) -> str:
        """Return a human description of the identified variant."""
        raise NotImplementedError
