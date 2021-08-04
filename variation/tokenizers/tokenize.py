"""A module for tokenizing."""
import re
from typing import Iterable, List
from .amplification import Amplification
from .deletion import Deletion
from .exon import Exon
from .expression import Expression
from .fusion import Fusion
from .gain_of_function import GainOfFunction
from .gene_pair import GenePair
from .gene_symbol import GeneSymbol
from .loss_of_function import LossOfFunction
from .overexpression import OverExpression
from .protein_alternate import ProteinAlternate
from .protein_frameshift import ProteinFrameshift
from .protein_termination import ProteinTermination
from .underexpression import UnderExpression
from .amino_acid_substitution import AminoAcidSubstitution
from .polypeptide_truncation import PolypeptideTruncation
from .silent_mutation import SilentMutation
from .coding_dna_substitution import CodingDNASubstitution
from .genomic_substitution import GenomicSubstitution
from .coding_dna_silent_mutation import CodingDNASilentMutation
from .genomic_silent_mutation import GenomicSilentMutation
from .amino_acid_delins import AminoAcidDelIns
from .coding_dna_delins import CodingDNADelIns
from .genomic_delins import GenomicDelIns
from .wild_type import WildType
from .hgvs import HGVS
from .reference_sequence import ReferenceSequence
from .locus_reference_genomic import LocusReferenceGenomic
from .amino_acid_deletion import AminoAcidDeletion
from .coding_dna_deletion import CodingDNADeletion
from .genomic_deletion import GenomicDeletion
from .amino_acid_insertion import AminoAcidInsertion
from .coding_dna_insertion import CodingDNAInsertion
from .genomic_insertion import GenomicInsertion
from variation.schemas.token_response_schema import Token, TokenMatchType
from .caches import GeneSymbolCache, AminoAcidCache, NucleotideCache
from variation import HGNC_GENE_SYMBOL_PATH


class Tokenize:
    """The tokenize class."""

    def __init__(self, gene_file_path=HGNC_GENE_SYMBOL_PATH) -> None:
        """Initialize the tokenize class.

        :param str gene_file_path: The path to the gene file
        """
        gene_cache = GeneSymbolCache(gene_file_path)
        amino_acid_cache = AminoAcidCache()
        nucleotide_cache = NucleotideCache()

        self.search_term_splitter = re.compile(r'\s+')

        self.tokenizers = (
            Amplification(),
            Deletion(),
            Exon(),
            Expression(),
            Fusion(),
            GainOfFunction(),
            GenePair(gene_cache),
            GeneSymbol(gene_cache),
            LossOfFunction(),
            OverExpression(),
            ProteinAlternate(amino_acid_cache),
            ProteinFrameshift(amino_acid_cache),
            AminoAcidSubstitution(amino_acid_cache),
            PolypeptideTruncation(amino_acid_cache),
            SilentMutation(amino_acid_cache),
            CodingDNASubstitution(),
            GenomicSubstitution(),
            CodingDNASilentMutation(),
            GenomicSilentMutation(),
            AminoAcidDelIns(amino_acid_cache, nucleotide_cache),
            CodingDNADelIns(amino_acid_cache, nucleotide_cache),
            GenomicDelIns(amino_acid_cache, nucleotide_cache),
            AminoAcidDeletion(amino_acid_cache, nucleotide_cache),
            CodingDNADeletion(amino_acid_cache, nucleotide_cache),
            GenomicDeletion(amino_acid_cache, nucleotide_cache),
            AminoAcidInsertion(amino_acid_cache, nucleotide_cache),
            CodingDNAInsertion(amino_acid_cache, nucleotide_cache),
            GenomicInsertion(amino_acid_cache, nucleotide_cache),
            ProteinTermination(amino_acid_cache),
            UnderExpression(),
            WildType(),
            HGVS(),
            ReferenceSequence(),
            LocusReferenceGenomic(),
        )

    def perform(self, search_string: str, warnings: List[str])\
            -> Iterable[Token]:
        """Return an iterable of tokens for a given search string.

        :param str search_string: The input string to search on
        :param list warnings: List of warnings
        :return: An Iterable of Tokens
        """
        tokens: List[Token] = list()
        terms = self.search_term_splitter.split(search_string)
        self._add_tokens(tokens, terms, search_string, warnings)

        # If reference sequence: Check description
        if list(map(lambda t: t.token_type, tokens)) == ['ReferenceSequence']:
            try:
                self._add_tokens(tokens,
                                 [search_string.split(':')[1]],
                                 search_string, warnings)
            except IndexError:
                return tokens
        return tokens

    def _add_tokens(self, tokens, terms, search_string, warnings):
        """Add tokens to a list for a given search string.

        :param list tokens: A list of tokens
        :param str search_string: The input string to search on
        :param list warnings: List of warnings
        """
        for term in terms:
            if not term:
                continue
            matched = False
            for tokenizer in self.tokenizers:
                res = tokenizer.match(term)
                if res:
                    tokens.append(res)
                    token = list(map(lambda t: t.token_type, tokens))[0]
                    if token == 'HGVS' or token == 'LocusReferenceGenomic':
                        # Give specific type of HGVS (i.e. protein sub)
                        if len(tokens) == 1:
                            self._add_tokens(tokens,
                                             [search_string.split(':')[1]],
                                             search_string, warnings)
                    matched = True
                    break
                else:
                    continue
            if not matched:
                warnings.append(f"Unable to tokenize {term}")
                tokens.append(Token(
                    token='',
                    token_type='Unknown',
                    input_string=term,
                    match_type=TokenMatchType.UNSPECIFIED
                ))
