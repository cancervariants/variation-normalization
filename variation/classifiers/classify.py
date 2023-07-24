"""Module for classification."""
from typing import List

from variation.schemas.classification_response_schema import Classification
from variation.schemas.token_response_schema import Token, TokenType
from variation.classifiers import (
    ProteinSubstitutionClassifier, ProteinReferenceAgreeClassifier,
    CdnaSubstitutionClassifier, GenomicSubstitutionClassifier,
    CdnaReferenceAgreeClassifier, GenomicReferenceAgreeClassifier,
    ProteinDelInsClassifier, CdnaDelInsClassifier, GenomicDelInsClassifier,
    ProteinDeletionClassifier, CdnaDeletionClassifier, GenomicDeletionClassifier,
    ProteinInsertionClassifier, CdnaInsertionClassifier,
    GenomicInsertionClassifier, GenomicDuplicationAmbiguousClassifier,
    GenomicDuplicationClassifier, GenomicDeletionAmbiguousClassifier,
    ProteinStopGainClassifier, AmplificationClassifier, HgvsClassifier,
    GnomadVcfClassifier, Classifier
)


class Classify:
    """The classify class."""

    hgvs_classifier = HgvsClassifier()
    gnomad_vcf_classifier = GnomadVcfClassifier()
    classifiers: List[Classifier] = [
        ProteinDelInsClassifier(),
        ProteinSubstitutionClassifier(),
        ProteinStopGainClassifier(),
        ProteinReferenceAgreeClassifier(),
        CdnaSubstitutionClassifier(),
        GenomicSubstitutionClassifier(),
        CdnaReferenceAgreeClassifier(),
        GenomicReferenceAgreeClassifier(),
        ProteinDelInsClassifier(),
        CdnaDelInsClassifier(),
        GenomicDelInsClassifier(),
        ProteinDeletionClassifier(),
        CdnaDeletionClassifier(),
        GenomicDeletionClassifier(),
        GenomicDeletionAmbiguousClassifier(),
        ProteinInsertionClassifier(),
        CdnaInsertionClassifier(),
        GenomicInsertionClassifier(),
        GenomicDuplicationClassifier(),
        GenomicDuplicationAmbiguousClassifier(),
        AmplificationClassifier()
    ]

    def perform(self, tokens: List[Token]) -> Classification:
        """Classify a list of tokens.

        :param tokens: List of tokens found
        :return: Classification for a list of tokens
        """
        classification = None

        if len(tokens) == 1:
            token_type = tokens[0].token_type

            if token_type == TokenType.HGVS:
                classification = self.hgvs_classifier.match(tokens[0])
            elif token_type == TokenType.GNOMAD_VCF:
                classification = self.gnomad_vcf_classifier.match(tokens[0])
        else:
            for classifier in self.classifiers:
                # We only do EXACT match candidates
                can_classify = classifier.can_classify(tokens)
                if can_classify:
                    classification = classifier.match(tokens)
                    if classification:
                        break

        return classification
