"""A module for the Cdna DelIns Classifier."""

from variation.classifiers.classifier import Classifier
from variation.schemas.classification_response_schema import (
    CdnaDelInsClassification,
    Nomenclature,
)
from variation.schemas.token_response_schema import Token, TokenType


class CdnaDelInsClassifier(Classifier):
    """The Cdna DelIns Classifier class."""

    def exact_match_candidates(self) -> list[list[TokenType]]:
        """Return the token match candidates for the cdna delins classification.

        :return: List of list of tokens, where order matters, that represent a cdna
        delins classification.
        """
        return [
            [TokenType.GENE, TokenType.CDNA_DELINS],
            [
                TokenType.GENE,
                TokenType.PROTEIN_DELETION,
                TokenType.CDNA_DELINS,
            ],
            [
                TokenType.GENE,
                TokenType.PROTEIN_DELINS,
                TokenType.CDNA_DELINS,
            ],
            [
                TokenType.GENE,
                TokenType.PROTEIN_FRAMESHIFT,
                TokenType.CDNA_DELINS,
            ],
            [
                TokenType.GENE,
                TokenType.PROTEIN_INSERTION,
                TokenType.CDNA_DELINS,
            ],
            [
                TokenType.GENE,
                TokenType.PROTEIN_REFERENCE_AGREE,
                TokenType.CDNA_DELINS,
            ],
            [
                TokenType.GENE,
                TokenType.PROTEIN_STOP_GAIN,
                TokenType.CDNA_DELINS,
            ],
            [
                TokenType.GENE,
                TokenType.PROTEIN_SUBSTITUTION,
                TokenType.CDNA_DELINS,
            ],
        ]

    def match(self, tokens: list[Token]) -> CdnaDelInsClassification:
        """Return the cdna delins classification from a list of token matches.

        :param tokens: List of ordered tokens that are exact match candidates for a
            cdna delins classification
        :return: cdna delins classification for the list of matched tokens
        """
        gene_token, cdna_delins_token = self.split_gene_cdna_tokens(tokens)

        return CdnaDelInsClassification(
            matching_tokens=tokens,
            nomenclature=Nomenclature.FREE_TEXT,
            gene_token=gene_token,
            pos0=cdna_delins_token.pos0,
            pos1=cdna_delins_token.pos1,
            inserted_sequence=cdna_delins_token.inserted_sequence,
        )
