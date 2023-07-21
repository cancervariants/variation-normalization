"""A module for the HGVS Classifier."""
from typing import Dict, List, Optional

from variation.schemas.classification_response_schema import (
    Classification, ClassificationType, GenomicSubstitutionClassification,
    GenomicInsertionClassification, ProteinSubstitutionClassification,
    CdnaSubstitutionClassification, CdnaDeletionClassification,
    ProteinDeletionClassification, ProteinDelInsClassification,
    ProteinInsertionClassification, ProteinReferenceAgreeClassification,
    CdnaDelInsClassification, GenomicDelInsClassification, CdnaInsertionClassification,
    CdnaReferenceAgreeClassification, GenomicReferenceAgreeClassification,
    ProteinStopGainClassification, GenomicDuplicationClassification,
    GenomicDuplicationAmbiguousClassification, GenomicDeletionClassification,
    GenomicDeletionAmbiguousClassification,
    Nomenclature, SequenceOntology
)
from variation.schemas.token_response_schema import HgvsToken, TokenType, CoordinateType
from variation.schemas.app_schemas import AmbiguousRegexType
from variation.regex import (
    GENOMIC_REGEXPRS, GENOMIC_DUP_AMBIGUOUS_REGEXPRS, PROTEIN_REGEXPRS, CDNA_REGEXPRS,
    GENOMIC_DEL_AMBIGUOUS_REGEXPRS
)
from variation.classifiers import Classifier
from variation.utils import get_ambiguous_type


class HgvsClassifier(Classifier):
    """The HGVS Classifier."""

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the exact match token type candidates."""
        return [
            [TokenType.HGVS]
        ]

    def match(self, token: HgvsToken) -> Optional[Classification]:
        classification = None
        params = {
            "matching_tokens": [token],
            "nomenclature": Nomenclature.HGVS,
            "ac": token.accession
        }

        if token.coordinate_type == CoordinateType.LINEAR_GENOMIC:
            classification = self._genomic_classification(token, params)
            if not classification:
                # Try ambiguous
                classification = self._genomic_ambiguous_classification(token, params)
        elif token.coordinate_type == CoordinateType.CODING_DNA:
            classification = self._cdna_classification(token, params)
        elif token.coordinate_type == CoordinateType.PROTEIN:
            classification = self._protein_classification(token, params)

        return classification

    def _protein_classification(
        self, token: HgvsToken, params: Dict
    ) -> Optional[Classification]:
        for regex, _, classification_type in PROTEIN_REGEXPRS:
            match = regex.match(token.change)

            if match:
                match_dict = match.groupdict()
                params.update(match_dict)

                if classification_type == ClassificationType.PROTEIN_SUBSTITUTION:
                    params["pos"] = int(params["pos"])
                    if params["alt"] in {"Ter", "*"}:
                        params["alt"] = "*"
                        return ProteinStopGainClassification(**params)
                    else:
                        return ProteinSubstitutionClassification(**params)
                elif classification_type == ClassificationType.PROTEIN_REFERENCE_AGREE:
                    params["pos"] = int(params["pos"])
                    return ProteinReferenceAgreeClassification(**params)
                elif classification_type == ClassificationType.PROTEIN_DELINS:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = int(params["pos1"]) if params["pos1"] is not None else params["pos1"]  # noqa: E501
                    return ProteinDelInsClassification(**params)
                elif classification_type == ClassificationType.PROTEIN_DELETION:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = int(params["pos1"]) if params["pos1"] is not None else params["pos1"]  # noqa: E501
                    return ProteinDeletionClassification(**params)
                elif classification_type == ClassificationType.PROTEIN_INSERTION:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = int(params["pos1"]) if params["pos1"] is not None else params["pos1"]  # noqa: E501
                    return ProteinInsertionClassification(**params)

        return None

    def _cdna_classification(
        self, token: HgvsToken, params: Dict
    ) -> Optional[Classification]:
        for regex, _, classification_type in CDNA_REGEXPRS:
            match = regex.match(token.change)

            if match:
                match_dict = match.groupdict()
                params.update(match_dict)

                if classification_type == ClassificationType.CODING_DNA_SUBSTITUTION:
                    len_ref = len(params["ref"])
                    len_alt = len(params["alt"])

                    if len_ref == 1 and len_alt == 1:
                        params["so_id"] = SequenceOntology.SNV
                    else:
                        params["so_id"] = SequenceOntology.MNV

                    params["pos"] = int(params["pos"])

                    return CdnaSubstitutionClassification(**params)
                elif classification_type == ClassificationType.CODING_DNA_REFERENCE_AGREE:  # noqa: E501
                    params["pos"] = int(params["pos"])
                    return CdnaReferenceAgreeClassification(**params)
                elif classification_type == ClassificationType.CODING_DNA_DELINS:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = int(params["pos1"]) if params["pos1"] is not None else params["pos1"]  # noqa: E501
                    return CdnaDelInsClassification(**params)
                elif classification_type == ClassificationType.CODING_DNA_DELETION:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = int(params["pos1"]) if params["pos1"] is not None else params["pos1"]  # noqa: E501
                    return CdnaDeletionClassification(**params)
                elif classification_type == ClassificationType.CODING_DNA_INSERTION:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = int(params["pos1"]) if params["pos1"] is not None else params["pos1"]  # noqa: E501
                    return CdnaInsertionClassification(**params)

    def _genomic_classification(
        self, token: HgvsToken, params: Dict
    ) -> Optional[Classification]:
        for regex, _, classification_type in GENOMIC_REGEXPRS:
            match = regex.match(token.change)

            if match:
                match_dict = match.groupdict()
                params.update(match_dict)

                if classification_type == ClassificationType.GENOMIC_SUBSTITUTION:
                    len_ref = len(params["ref"])
                    len_alt = len(params["alt"])

                    if len_ref == len_alt == 1:
                        params["so_id"] = SequenceOntology.SNV
                    else:
                        params["so_id"] = SequenceOntology.MNV

                    params["pos"] = int(params["pos"])

                    return GenomicSubstitutionClassification(**params)
                elif classification_type == ClassificationType.GENOMIC_REFERENCE_AGREE:
                    params["pos"] = int(params["pos"])
                    return GenomicReferenceAgreeClassification(**params)
                elif classification_type == ClassificationType.GENOMIC_DELINS:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = int(params["pos1"]) if params["pos1"] is not None else params["pos1"]  # noqa: E501
                    return GenomicDelInsClassification(**params)
                elif classification_type == ClassificationType.GENOMIC_INSERTION:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = int(params["pos1"]) if params["pos1"] is not None else params["pos1"]  # noqa: E501
                    return GenomicInsertionClassification(**params)
                elif classification_type == ClassificationType.GENOMIC_DELETION:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = int(params["pos1"]) if params["pos1"] is not None else params["pos1"]  # noqa: E501
                    return GenomicDeletionClassification(**params)
                elif classification_type == ClassificationType.GENOMIC_DUPLICATION:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = int(params["pos1"]) if params["pos1"] is not None else params["pos1"]  # noqa: E501
                    return GenomicDuplicationClassification(**params)

    def _genomic_ambiguous_classification(
        self, token: HgvsToken, params: Dict
    ) -> Optional[Classification]:
        if token.token.endswith("dup"):
            return self._genomic_dup_ambiguous_classification(token, params)
        elif token.token.endswith("del"):
            return self._genomic_del_ambiguous_classification(token, params)

        return None

    @staticmethod
    def _update_ambiguous_params(
        params: Dict, regex_type: AmbiguousRegexType
    ) -> None:
        params["pos0"] = int(params["pos0"]) if params["pos0"] != "?" else params["pos0"]  # noqa: E501

        if "pos1" in params:
            params["pos1"] = int(params["pos1"]) if params["pos1"] != "?" else params["pos1"]  # noqa: E501
        else:
            params["pos1"] = None

        params["pos2"] = int(params["pos2"]) if params["pos2"] != "?" else params["pos2"]  # noqa: E501

        if "pos3" in params:
            params["pos3"] = int(params["pos3"]) if params["pos3"] != "?" else params["pos3"]  # noqa: E501
        else:
            params["pos3"] = None

        ambiguous_type = get_ambiguous_type(
            params["pos0"], params["pos1"], params["pos2"], params["pos3"], regex_type
        )
        if ambiguous_type:
            params["ambiguous_type"] = ambiguous_type

    def _genomic_dup_ambiguous_classification(
        self, token: HgvsToken, params: Dict
    ) -> Optional[Classification]:
        for regex, _, classification_type, regex_type in GENOMIC_DUP_AMBIGUOUS_REGEXPRS:
            match = regex.match(token.change)

            if match:
                match_dict = match.groupdict()
                params.update(match_dict)

                if classification_type == ClassificationType.GENOMIC_DUPLICATION_AMBIGUOUS:  # noqa: E501
                    self._update_ambiguous_params(params, regex_type)

                    # If ambiguous type not in params, it means we don't support it yet
                    if "ambiguous_type" in params:
                        return GenomicDuplicationAmbiguousClassification(**params)
        return None

    def _genomic_del_ambiguous_classification(
        self, token: HgvsToken, params: Dict
    ) -> Optional[Classification]:
        for regex, _, classification_type, regex_type in GENOMIC_DEL_AMBIGUOUS_REGEXPRS:
            match = regex.match(token.change)

            if match:
                match_dict = match.groupdict()
                params.update(match_dict)

                if classification_type == ClassificationType.GENOMIC_DELETION_AMBIGUOUS:  # noqa: E501
                    self._update_ambiguous_params(params, regex_type)

                    # If ambiguous type not in params, it means we don't support it yet
                    if "ambiguous_type" in params:
                        return GenomicDeletionAmbiguousClassification(**params)
        return None
