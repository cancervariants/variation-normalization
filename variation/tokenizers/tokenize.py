"""A module for tokenizing."""
from typing import Iterable, List

from variation.schemas.token_response_schema import Token, TokenType
from .gene_symbol import GeneSymbol
from .protein_substitution import ProteinSubstitution
from .protein_reference_agree import ProteinReferenceAgree
from .cdna_and_genomic_reference_agree import CdnaGenomicReferenceAgree
from .cdna_substitution import CdnaSubstitution
from .genomic_substitution import GenomicSubstitution
from .protein_delins import ProteinDelIns
from .cdna_delins import CdnaDelIns
from .genomic_delins import GenomicDelIns
from .hgvs import HGVS
from .protein_deletion import ProteinDeletion
from .cdna_deletion import CdnaDeletion
from .genomic_deletion import GenomicDeletion
from .protein_insertion import ProteinInsertion
from .cdna_insertion import CdnaInsertion
from .genomic_insertion import GenomicInsertion
from .genomic_duplication import GenomicDuplication
from .gnomad_vcf import GnomadVCF
from .free_text_categorical import FreeTextCategorical

r"(\((\?|d+)_(\?|\d+)\))_(\((\?|\d+)_(\?|\d+)\))dup"


class Tokenize:
    """The tokenize class."""

    def __init__(self, gene_symbol: GeneSymbol) -> None:
        """Initialize the tokenize class."""
        self.gene_symbol = gene_symbol
        self.tokenizers = [
            HGVS(),
            GnomadVCF(),
            self.gene_symbol,
            FreeTextCategorical(),
            # Substitution
            ProteinSubstitution(),
            GenomicSubstitution(),
            CdnaSubstitution(),
            # Reference Agree
            ProteinReferenceAgree(),
            CdnaGenomicReferenceAgree(),
            # Delins
            ProteinDelIns(),
            CdnaDelIns(),
            GenomicDelIns(),
            # Deletion
            ProteinDeletion(),
            CdnaDeletion(),
            GenomicDeletion(),
            # Insertion
            ProteinInsertion(),
            CdnaInsertion(),
            GenomicInsertion(),
            # Duplication
            GenomicDuplication()
        ]

    def perform(self, search_string: str, warnings: List[str]) -> Iterable[Token]:
        """Return an iterable of tokens for a given search string.

        :param str search_string: The input string to search on
        :param List warnings: List of warnings
        :return: An Iterable of Tokens
        """
        terms = search_string.split()

        tokens: List[Token] = []
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
                        matched = True
                        break

            if not matched:
                warnings.append(f"Unable to tokenize: {term}")
                tokens.append(Token(
                    token=term,
                    token_type=TokenType.UNKNOWN,
                    input_string=term
                ))

        return tokens
