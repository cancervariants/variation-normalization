"""A module for the gnomAD VCF Classifier"""
from typing import List, Optional, Union

from variation.classifiers.classifier import Classifier
from variation.schemas.classification_response_schema import (
    GenomicDeletionClassification,
    GenomicDelInsClassification,
    GenomicInsertionClassification,
    GenomicReferenceAgreeClassification,
    GenomicSubstitutionClassification,
    Nomenclature,
)
from variation.schemas.token_response_schema import GnomadVcfToken, TokenType


class GnomadVcfClassifier(Classifier):
    """The gnomAD VCF Classifier"""

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the token match candidates for the gnomad vcf classification.

        :return: List of list of tokens, where order matters, that represent a gnomad
        vcf classification.
        """
        return [[TokenType.GNOMAD_VCF]]

    def match(
        self, token: GnomadVcfToken
    ) -> Optional[
        Union[
            GenomicReferenceAgreeClassification,
            GenomicSubstitutionClassification,
            GenomicInsertionClassification,
            GenomicDeletionClassification,
        ]
    ]:
        """Return the genomic classification (either reference agree, substitution,
        insertion, or deletion) from a gnomad vcf token.
        Currently only support simple genomic variation.

        :param token: gnomad vcf token
        :return: The corresponding genomic classification for the gnomad vcf token if
            simple variation change. Else, return `None`
        """
        params = {"matching_tokens": [token], "nomenclature": Nomenclature.GNOMAD_VCF}

        ref = token.ref
        alt = token.alt

        len_ref = len(ref)
        len_alt = len(alt)

        if len_ref == len_alt:
            # substitution
            params["pos"] = token.pos

            if ref == alt:
                return GenomicReferenceAgreeClassification(**params)
            else:
                params["ref"] = ref
                params["alt"] = alt

                return GenomicSubstitutionClassification(**params)

        # delins
        params["pos0"] = token.pos
        params["pos1"] = (params["pos0"] + len_ref) - 1
        if params["pos0"] == params["pos1"]:
            del params["pos1"]

        params["inserted_sequence"] = alt
        return GenomicDelInsClassification(**params)
