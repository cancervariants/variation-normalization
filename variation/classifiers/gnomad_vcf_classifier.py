"""A module for the gnomAD VCF Classifier"""
import re
from typing import List, Optional, Union

from variation.classifiers.classifier import Classifier
from variation.schemas.classification_response_schema import (
    GenomicDeletionClassification,
    GenomicDelInsClassification,
    GenomicInsertionClassification,
    GenomicReferenceAgreeClassification,
    GenomicSubstitutionClassification,
    Nomenclature,
    SequenceOntology,
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

    @staticmethod
    def _get_substring_end_index(pattern: str, target_str: str) -> int:
        """Get end index of the pattern if the pattern starts at the beginning of
        the string to match

        :param pattern: Pattern to find at the beginning of `string`
        :param target_str: The target string
        :return: End index of match if `pattern` found at the beginning of `target_str`.
            If no match found or match is not found at the beginning of `target_str`,
            returns 0
        """
        matches = [m for m in re.finditer(pattern, target_str)]
        match_end = 0
        if matches:
            match = matches[0]
            if match.start() == 0:
                match_end = match.end()
        return match_end

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
                params["so_id"] = SequenceOntology.NO_SEQUENCE_ALTERATION
                return GenomicReferenceAgreeClassification(**params)
            else:
                if len_ref == 1:
                    params["so_id"] = SequenceOntology.SNV
                else:
                    params["so_id"] = SequenceOntology.MNV

                params["ref"] = ref
                params["alt"] = alt

                return GenomicSubstitutionClassification(**params)
        elif len_ref < len_alt:
            if len_ref == 1 and alt[0] == ref:
                params["pos0"] = token.pos
                params["pos1"] = params["pos0"] + 1
                params["inserted_sequence"] = alt[1:]
                return GenomicInsertionClassification(**params)
        else:
            match_end = self._get_substring_end_index(alt, ref)
            if match_end:
                params["pos0"] = token.pos + match_end
                params["deleted_sequence"] = ref[match_end:]
                params["pos1"] = params["pos0"] + len(params["deleted_sequence"])
                return GenomicDeletionClassification(**params)

        # delins
        params["pos0"] = token.pos
        params["pos1"] = (params["pos0"] + len_ref) - 1
        if params["pos0"] == params["pos1"]:
            del params["pos1"]

        params["inserted_sequence"] = alt
        return GenomicDelInsClassification(**params)
