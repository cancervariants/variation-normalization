"""A module for tokenizing."""
from typing import Iterable, List

from variation.schemas.token_response_schema import Token, TokenMatchType
from variation.tokenizers.caches.amino_acid_cache import AminoAcidCache
from .gene_symbol import GeneSymbol
from .protein_substitution import ProteinSubstitution
from .polypeptide_truncation import PolypeptideTruncation
from .silent_mutation import SilentMutation
from .coding_dna_substitution import CodingDNASubstitution
from .genomic_substitution import GenomicSubstitution
from .coding_dna_silent_mutation import CodingDNASilentMutation
from .genomic_silent_mutation import GenomicSilentMutation
from .protein_delins import ProteinDelIns
from .coding_dna_delins import CodingDNADelIns
from .genomic_delins import GenomicDelIns
from .hgvs import HGVS
from .reference_sequence import ReferenceSequence
from .locus_reference_genomic import LocusReferenceGenomic
from .protein_deletion import ProteinDeletion
from .coding_dna_deletion import CodingDNADeletion
from .genomic_deletion import GenomicDeletion
from .protein_insertion import ProteinInsertion
from .coding_dna_insertion import CodingDNAInsertion
from .genomic_insertion import GenomicInsertion
from .genomic_uncertain_deletion import GenomicUncertainDeletion
from .genomic_duplication import GenomicDuplication
from .genomic_deletion_range import GenomicDeletionRange
from .gnomad_vcf import GnomadVCF
from .caches import NucleotideCache


class Tokenize:
    """The tokenize class."""

    def __init__(self, amino_acid_cache: AminoAcidCache,
                 gene_symbol: GeneSymbol) -> None:
        """Initialize the tokenize class."""
        nucleotide_cache = NucleotideCache()
        self.tokenizers = (
            HGVS(),
            ReferenceSequence(),
            LocusReferenceGenomic(),
            GnomadVCF(),
            gene_symbol,
            ProteinSubstitution(amino_acid_cache),
            PolypeptideTruncation(amino_acid_cache),
            SilentMutation(amino_acid_cache),
            CodingDNASubstitution(),
            GenomicSubstitution(),
            CodingDNASilentMutation(),
            GenomicSilentMutation(),
            ProteinDelIns(amino_acid_cache, nucleotide_cache),
            CodingDNADelIns(amino_acid_cache, nucleotide_cache),
            GenomicDelIns(amino_acid_cache, nucleotide_cache),
            ProteinDeletion(amino_acid_cache, nucleotide_cache),
            CodingDNADeletion(amino_acid_cache, nucleotide_cache),
            GenomicDeletion(amino_acid_cache, nucleotide_cache),
            ProteinInsertion(amino_acid_cache, nucleotide_cache),
            CodingDNAInsertion(amino_acid_cache, nucleotide_cache),
            GenomicInsertion(amino_acid_cache, nucleotide_cache),
            GenomicUncertainDeletion(),
            GenomicDuplication(),
            GenomicDeletionRange()
        )

    def perform(self, search_string: str, warnings: List[str]) -> Iterable[Token]:
        """Return an iterable of tokens for a given search string.

        :param str search_string: The input string to search on
        :param List warnings: List of warnings
        :return: An Iterable of Tokens
        """
        tokens: List[Token] = list()
        # Currently splits on whitespace
        # Also adds a token if an string looks like an accession
        #  ex: NC_, NM_, ENST
        terms = search_string.split()
        self._add_tokens(tokens, terms, search_string, warnings)
        return tokens

    def _add_tokens(self, tokens: List, terms: List, search_string: str,
                    warnings: List) -> None:
        """Add tokens to a list for a given search string.

        :param List tokens: A list of tokens
        :param List terms: List of terms from query
        :param str search_string: The input string to search on
        :param List warnings: List of warnings
        """
        for term in terms:
            if not term:
                continue
            matched = False
            for tokenizer in self.tokenizers:
                res = tokenizer.match(term)
                if res:
                    if isinstance(res, List):
                        for r in res:
                            tokens.append(r)
                            if not matched:
                                matched = True
                    else:
                        tokens.append(res)
                        token = list(map(lambda t: t.token_type, tokens))[0]
                        if token == "HGVS" or \
                                token == "LocusReferenceGenomic" \
                                or token == "ReferenceSequence":
                            # Give specific type of HGVS (i.e. protein sub)
                            if len(tokens) == 1:
                                self._add_tokens(tokens,
                                                 [search_string.split(":")[1]],
                                                 search_string, warnings)
                        matched = True
                        break
                else:
                    continue
            if not matched:
                warnings.append(f"Unable to tokenize {term}")
                tokens.append(Token(
                    token="",
                    token_type="Unknown",
                    input_string=term,
                    match_type=TokenMatchType.UNSPECIFIED
                ))
