"""Module for classification."""
from typing import List, Callable
from variation.schemas.classification_response_schema import Classification, \
    ConfidenceRating
from variation.schemas.token_response_schema import Token
from variation.classifiers import ComplexClassifier, ExpressionClassifier, \
    FusionClassifier, OncogenicClassifier, ProteinAlternateClassifier, \
    ProteinDelinsClassifier, ProteinFrameshiftClassifier, \
    ProteinTerminationClassifier, \
    AminoAcidSubstitutionClassifier, PolypeptideTruncationClassifier, \
    SilentMutationClassifier, CodingDNASubstitutionClassifier,\
    GenomicSubstitutionClassifier, CodingDNASilentMutationClassifier, \
    GenomicSilentMutationClassifier, AminoAcidDelInsClassifier, \
    CodingDNADelInsClassifier, GenomicDelInsClassifier, \
    AminoAcidDeletionClassifier, CodingDNADeletionClassifier, \
    GenomicDeletionClassifier, AminoAcidInsertionClassifier, \
    CodingDNAInsertionClassifier, GenomicInsertionClassifier, Classifier


class Classify:
    """The classify class."""

    def __init__(self) -> None:
        """Initialize the Classify class."""
        self.classifiers: List[Classifier] = [
            ComplexClassifier(),
            ExpressionClassifier(),
            FusionClassifier(),
            OncogenicClassifier(),
            ProteinAlternateClassifier(),
            ProteinDelinsClassifier(),
            ProteinFrameshiftClassifier(),
            AminoAcidSubstitutionClassifier(),
            PolypeptideTruncationClassifier(),
            SilentMutationClassifier(),
            ProteinTerminationClassifier(),
            CodingDNASubstitutionClassifier(),
            GenomicSubstitutionClassifier(),
            CodingDNASilentMutationClassifier(),
            GenomicSilentMutationClassifier(),
            AminoAcidDelInsClassifier(),
            CodingDNADelInsClassifier(),
            GenomicDelInsClassifier(),
            AminoAcidDeletionClassifier(),
            CodingDNADeletionClassifier(),
            GenomicDeletionClassifier(),
            AminoAcidInsertionClassifier(),
            CodingDNAInsertionClassifier(),
            GenomicInsertionClassifier()
        ]

    def perform(self, tokens: List[Token]) -> List[Classification]:
        """Classify a list of tokens."""
        classifications: List[Classification] = list()

        for classifier in self.classifiers:
            res = classifier.match(tokens)
            if res is not None and res.confidence == ConfidenceRating.EXACT:
                return [res]
            elif res is not None:
                classifications.append(res)

        highest_confidence = 0
        for match in classifications:
            if match.confidence.value > highest_confidence:
                highest_confidence = match.confidence.value

        filter_lambda: Callable[[Classification], bool] = \
            lambda match: match.confidence.value == highest_confidence

        return list(filter(filter_lambda, classifications))
