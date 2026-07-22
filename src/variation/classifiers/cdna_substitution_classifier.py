"""A module for the Cdna Substitution Classifier."""

from variation.classifiers.classifier import Classifier
from variation.schemas.classification_response_schema import (
    CdnaSubstitutionClassification,
    Nomenclature,
)
from variation.schemas.token_response_schema import Token, TokenType


class CdnaSubstitutionClassifier(Classifier):
    """The Cdna Substitution Classifier class."""

    def exact_match_candidates(self) -> list[list[TokenType]]:
        """Return the token match candidates for the cdna substitution classification.

        :return: List of list of tokens, where order matters, that represent a cdna
        substitution classification.
        """
        return [
            [TokenType.GENE, TokenType.CDNA_SUBSTITUTION],
            [
                TokenType.GENE,
                TokenType.PROTEIN_DELETION,
                TokenType.CDNA_SUBSTITUTION,
            ],
            [
                TokenType.GENE,
                TokenType.PROTEIN_DELINS,
                TokenType.CDNA_SUBSTITUTION,
            ],
            [
                TokenType.GENE,
                TokenType.PROTEIN_FRAMESHIFT,
                TokenType.CDNA_SUBSTITUTION,
            ],
            [
                TokenType.GENE,
                TokenType.PROTEIN_INSERTION,
                TokenType.CDNA_SUBSTITUTION,
            ],
            [
                TokenType.GENE,
                TokenType.PROTEIN_REFERENCE_AGREE,
                TokenType.CDNA_SUBSTITUTION,
            ],
            [
                TokenType.GENE,
                TokenType.PROTEIN_STOP_GAIN,
                TokenType.CDNA_SUBSTITUTION,
            ],
            [
                TokenType.GENE,
                TokenType.PROTEIN_SUBSTITUTION,
                TokenType.CDNA_SUBSTITUTION,
            ],
        ]

    def match(self, tokens: list[Token]) -> CdnaSubstitutionClassification | None:
        """Return the cdna substitution classification from a list of token matches.

        :param tokens: List of ordered tokens that are exact match candidates for a
            cdna substitution classification
        :return: cdna substitution classification for the list of matched tokens
        """
        gene_token, cdna_sub_token = self.split_gene_cdna_tokens(tokens)

        if len(cdna_sub_token.ref) == len(cdna_sub_token.alt):
            return CdnaSubstitutionClassification(
                matching_tokens=tokens,
                nomenclature=Nomenclature.FREE_TEXT,
                gene_token=gene_token,
                pos=cdna_sub_token.pos,
                ref=cdna_sub_token.ref,
                alt=cdna_sub_token.alt,
            )

        return None
