"""Module for classification."""
from typing import List, Callable

from variation.schemas.classification_response_schema import Classification
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import (
    ProteinSubstitutionClassifier,
    ProteinReferenceAgreeClassifier, CodingDNASubstitutionClassifier,
    GenomicSubstitutionClassifier, CodingDNAReferenceAgreeClassifier,
    GenomicReferenceAgreeClassifier, ProteinDelInsClassifier, CodingDNADelInsClassifier,
    GenomicDelInsClassifier, ProteinDeletionClassifier, CdnaDeletionClassifier,
    GenomicDeletionClassifier, ProteinInsertionClassifier, CodingDNAInsertionClassifier,
    GenomicInsertionClassifier, GenomicUncertainDeletionClassifier,
    GenomicDuplicationClassifier, GenomicDeletionRangeClassifier,
    AmplificationClassifier, HgvsClassifier, Classifier
)


class Classify:
    """The classify class."""

    hgvs_classifier = HgvsClassifier()

    def __init__(self) -> None:
        """Initialize the Classify class."""
        self.classifiers: List[Classifier] = [
            # ProteinDelInsClassifier(),
            # ProteinSubstitutionClassifier(),
            # ProteinReferenceAgreeClassifier(),
            CodingDNASubstitutionClassifier(),
            # GenomicSubstitutionClassifier(),
            # CodingDNAReferenceAgreeClassifier(),
            # GenomicReferenceAgreeClassifier(),
            # ProteinDelInsClassifier(),
            # CodingDNADelInsClassifier(),
            # GenomicDelInsClassifier(),
            # ProteinDeletionClassifier(),
            CdnaDeletionClassifier(),
            # GenomicDeletionClassifier(),
            # ProteinInsertionClassifier(),
            # CodingDNAInsertionClassifier(),
            # GenomicInsertionClassifier(),
            # GenomicUncertainDeletionClassifier(),
            # GenomicDuplicationClassifier(),
            # GenomicDeletionRangeClassifier(),
            AmplificationClassifier()
        ]

    def perform(self, tokens: List[Token]) -> List[Classification]:
        """Classify a list of tokens.

        :param tokens: List of tokens found
        :return: Classifications
        """
        classifications: List[Classification] = list()

        if len(tokens) == 1 and tokens[0].token_type == TokenType.HGVS:
            res = self.hgvs_classifier.match(tokens[0])
            if res:
                classifications.append(res)
        else:
            # TODO: Still need to do this
            for classifier in self.classifiers:
                # We only do EXACT match candidates
                can_classify = classifier.can_classify(tokens)
                if can_classify:
                    res = classifier.match(tokens)
                    if res:
                        classifications.append(res)

        return classifications
