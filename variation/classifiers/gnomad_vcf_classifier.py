"""A module for the gnomAD VCF Classifier"""
from typing import List, Optional

from variation.schemas.classification_response_schema import (
    Classification, GenomicSubstitutionClassification, Nomenclature, SequenceOntology,
    GenomicReferenceAgreeClassification, GenomicInsertionClassification,
    GenomicDeletionClassification, GenomicDelInsClassification
)
from variation.schemas.token_response_schema import (
    GnomadVcfToken, TokenType
)
from variation.classifiers import Classifier


class GnomadVcfClassifier(Classifier):
    """The gnomAD VCF Classifier"""

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.GNOMAD_VCF]
        ]

    def match(self, token: GnomadVcfToken) -> Optional[Classification]:
        params = {
            "matching_tokens": [token],
            "nomenclature": Nomenclature.GNOMAD_VCF
        }

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
            # insertion
            if len_ref == 1:
                if ref[0] == alt[0]:
                    params["pos0"] = token.pos
                    params["pos1"] = token.pos + 1
                    params["inserted_sequence"] = alt[1:]
                    return GenomicInsertionClassification(**params)
        else:
            # deletion
            if len_alt == 1:
                if ref[0] == alt[0]:
                    del_seq = ref[1:]
                    params["pos0"] = token.pos
                    params["pos1"] = params["pos0"] + len(del_seq)
                    params["deleted_sequence"] = del_seq
                    return GenomicDeletionClassification(**params)

        return None
