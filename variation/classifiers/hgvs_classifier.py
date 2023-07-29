"""A module for the HGVS Classifier."""
from typing import Dict, List, Optional
from re import Match, Pattern

from variation.schemas.classification_response_schema import (
    Classification,
    ClassificationType,
    GenomicSubstitutionClassification,
    GenomicInsertionClassification,
    ProteinSubstitutionClassification,
    CdnaSubstitutionClassification,
    CdnaDeletionClassification,
    ProteinDeletionClassification,
    ProteinDelInsClassification,
    ProteinInsertionClassification,
    ProteinReferenceAgreeClassification,
    CdnaDelInsClassification,
    GenomicDelInsClassification,
    CdnaInsertionClassification,
    CdnaReferenceAgreeClassification,
    GenomicReferenceAgreeClassification,
    ProteinStopGainClassification,
    GenomicDuplicationClassification,
    GenomicDuplicationAmbiguousClassification,
    GenomicDeletionClassification,
    GenomicDeletionAmbiguousClassification,
    Nomenclature,
    SequenceOntology,
)
from variation.schemas.token_response_schema import HgvsToken, TokenType, CoordinateType
from variation.schemas.app_schemas import AmbiguousRegexType
from variation.regex import (
    GENOMIC_REGEXPRS,
    GENOMIC_DUP_AMBIGUOUS_REGEXPRS,
    PROTEIN_REGEXPRS,
    CDNA_REGEXPRS,
    GENOMIC_DEL_AMBIGUOUS_REGEXPRS,
)
from variation.classifiers import Classifier
from variation.utils import get_ambiguous_type


class HgvsClassifier(Classifier):
    """The HGVS Classifier."""

    def exact_match_candidates(self) -> List[List[TokenType]]:
        """Return the token match candidates for the hgvs classification.

        :return: List of list of tokens, where order matters, that represent a hgvs
        classification.
        """
        return [[TokenType.HGVS]]

    def match(self, token: HgvsToken) -> Optional[Classification]:
        """Return the classification from a hgvs token using regex matches to determine
        the type of classification.

        :param token: hgvs token
        :return: The corresponding classification for the hgvs token if a regex match
            is found. Else, return `None`
        """
        classification = None
        params = {
            "matching_tokens": [token],
            "nomenclature": Nomenclature.HGVS,
            "ac": token.accession,
        }

        if token.coordinate_type == CoordinateType.LINEAR_GENOMIC:
            classification = self._genomic_classification(token, params)
            if not classification:
                # Try ambiguous
                classification = self._genomic_ambiguous_classification(token, params)
        elif token.coordinate_type == CoordinateType.CDNA:
            classification = self._cdna_classification(token, params)
        elif token.coordinate_type == CoordinateType.PROTEIN:
            classification = self._protein_classification(token, params)

        return classification

    @staticmethod
    def _regex_match(change: str, regex: Pattern) -> Optional[Match]:
        """Strip parentheses from `change` and return whether or not `change` matches
        the `regex`

        :param change: The alteration part of the hgvs expression
        :param regex: The pattern to match against
        :return: A regex match if found against pattern, else `None`
        """
        if change[0] == "(" and change[-1] == ")":
            match = regex.match(change[1:-1])
        else:
            match = regex.match(change)
        return match

    def _protein_classification(
        self, token: HgvsToken, params: Dict
    ) -> Optional[Classification]:
        """Determine if hgvs token matches regex checks and return corresponding protein
        classification if a match is found

        :param token: hgvs token
        :param params: Base fields for a classification. This will get mutated if a
            match is found.
        :return: Protein classification if hgvs token matches regex checks. Else, `None`
        """
        for regex, _, classification_type in PROTEIN_REGEXPRS:
            match = self._regex_match(token.change, regex)

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
                    params["pos1"] = (
                        int(params["pos1"])
                        if params["pos1"] is not None
                        else params["pos1"]
                    )  # noqa: E501
                    return ProteinDelInsClassification(**params)
                elif classification_type == ClassificationType.PROTEIN_DELETION:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = (
                        int(params["pos1"])
                        if params["pos1"] is not None
                        else params["pos1"]
                    )  # noqa: E501
                    return ProteinDeletionClassification(**params)
                elif classification_type == ClassificationType.PROTEIN_INSERTION:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = (
                        int(params["pos1"])
                        if params["pos1"] is not None
                        else params["pos1"]
                    )  # noqa: E501
                    return ProteinInsertionClassification(**params)

        return None

    def _cdna_classification(
        self, token: HgvsToken, params: Dict
    ) -> Optional[Classification]:
        """Determine if hgvs token matches regex checks and return corresponding cdna
        classification if a match is found

        :param token: hgvs token
        :param params: Base fields for a classification. This will get mutated if a
            match is found.
        :return: cdna classification if hgvs token matches regex checks. Else, `None`
        """
        for regex, _, classification_type in CDNA_REGEXPRS:
            match = self._regex_match(token.change, regex)

            if match:
                match_dict = match.groupdict()
                params.update(match_dict)

                if classification_type == ClassificationType.CDNA_SUBSTITUTION:
                    len_ref = len(params["ref"])
                    len_alt = len(params["alt"])

                    if len_ref == 1 and len_alt == 1:
                        params["so_id"] = SequenceOntology.SNV
                    else:
                        params["so_id"] = SequenceOntology.MNV

                    params["pos"] = int(params["pos"])

                    return CdnaSubstitutionClassification(**params)
                elif (
                    classification_type == ClassificationType.CDNA_REFERENCE_AGREE
                ):  # noqa: E501
                    params["pos"] = int(params["pos"])
                    return CdnaReferenceAgreeClassification(**params)
                elif classification_type == ClassificationType.CDNA_DELINS:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = (
                        int(params["pos1"])
                        if params["pos1"] is not None
                        else params["pos1"]
                    )  # noqa: E501
                    return CdnaDelInsClassification(**params)
                elif classification_type == ClassificationType.CDNA_DELETION:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = (
                        int(params["pos1"])
                        if params["pos1"] is not None
                        else params["pos1"]
                    )  # noqa: E501
                    return CdnaDeletionClassification(**params)
                elif classification_type == ClassificationType.CDNA_INSERTION:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = (
                        int(params["pos1"])
                        if params["pos1"] is not None
                        else params["pos1"]
                    )  # noqa: E501
                    return CdnaInsertionClassification(**params)

    def _genomic_classification(
        self, token: HgvsToken, params: Dict
    ) -> Optional[Classification]:
        """Determine if hgvs token matches regex checks and return corresponding genomic
        classification if a match is found. Only checks against 'simple'
        duplication/deletions.

        :param token: hgvs token
        :param params: Base fields for a classification. This will get mutated if a
            match is found.
        :return: genomic classification if hgvs token matches regex checks. Else, `None`
        """
        for regex, _, classification_type in GENOMIC_REGEXPRS:
            match = self._regex_match(token.change, regex)

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
                    params["pos1"] = (
                        int(params["pos1"])
                        if params["pos1"] is not None
                        else params["pos1"]
                    )  # noqa: E501
                    return GenomicDelInsClassification(**params)
                elif classification_type == ClassificationType.GENOMIC_INSERTION:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = (
                        int(params["pos1"])
                        if params["pos1"] is not None
                        else params["pos1"]
                    )  # noqa: E501
                    return GenomicInsertionClassification(**params)
                elif classification_type == ClassificationType.GENOMIC_DELETION:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = (
                        int(params["pos1"])
                        if params["pos1"] is not None
                        else params["pos1"]
                    )  # noqa: E501
                    return GenomicDeletionClassification(**params)
                elif classification_type == ClassificationType.GENOMIC_DUPLICATION:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = (
                        int(params["pos1"])
                        if params["pos1"] is not None
                        else params["pos1"]
                    )  # noqa: E501
                    return GenomicDuplicationClassification(**params)

    def _genomic_ambiguous_classification(
        self, token: HgvsToken, params: Dict
    ) -> Optional[Classification]:
        """Determine if hgvs token matches regex checks and return corresponding genomic
        ambiguous classification if a match is found. Only checks against ambiguous
        duplication/deletions.

        :param token: hgvs token
        :param params: Base fields for a classification. This will get mutated if a
            match is found.
        :return: genomic ambiguous classification if hgvs token matches regex checks.
            Else, `None`
        """
        if token.token.endswith("dup"):
            return self._genomic_dup_ambiguous_classification(token, params)
        elif token.token.endswith("del"):
            return self._genomic_del_ambiguous_classification(token, params)

        return None

    @staticmethod
    def _update_ambiguous_params(params: Dict, regex_type: AmbiguousRegexType) -> None:
        """Mutates `params` to match correct types and gets associated ambiguous type
        from fields in `params`

        :param params: Fields for a classification. This will get mutated.
        :param regex_type: The kind of ambiguous regex that was used
        """
        params["pos0"] = (
            int(params["pos0"]) if params["pos0"] != "?" else params["pos0"]
        )  # noqa: E501

        if "pos1" in params:
            params["pos1"] = (
                int(params["pos1"]) if params["pos1"] != "?" else params["pos1"]
            )  # noqa: E501
        else:
            params["pos1"] = None

        params["pos2"] = (
            int(params["pos2"]) if params["pos2"] != "?" else params["pos2"]
        )  # noqa: E501

        if "pos3" in params:
            params["pos3"] = (
                int(params["pos3"]) if params["pos3"] != "?" else params["pos3"]
            )  # noqa: E501
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
        """Determine if hgvs token matches regex checks and return corresponding genomic
        ambiguous duplication classification if a match is found. Only checks against
        genomic ambiguous duplications.

        :param token: hgvs token
        :param params: Base fields for a classification. This will get mutated if a
            match is found.
        :return: genomic ambiguous duplication classification if hgvs token matches
            regex checks. Else, `None`
        """
        for regex, _, classification_type, regex_type in GENOMIC_DUP_AMBIGUOUS_REGEXPRS:
            match = regex.match(token.change)

            if match:
                match_dict = match.groupdict()
                params.update(match_dict)

                if (
                    classification_type
                    == ClassificationType.GENOMIC_DUPLICATION_AMBIGUOUS
                ):  # noqa: E501
                    self._update_ambiguous_params(params, regex_type)

                    # If ambiguous type not in params, it means we don't support it yet
                    if "ambiguous_type" in params:
                        return GenomicDuplicationAmbiguousClassification(**params)
        return None

    def _genomic_del_ambiguous_classification(
        self, token: HgvsToken, params: Dict
    ) -> Optional[Classification]:
        """Determine if hgvs token matches regex checks and return corresponding genomic
        ambiguous deletion classification if a match is found. Only checks against
        genomic ambiguous deletion.

        :param token: hgvs token
        :param params: Base fields for a classification. This will get mutated if a
            match is found.
        :return: genomic ambiguous deletion classification if hgvs token matches regex
            checks. Else, `None`
        """
        for regex, _, classification_type, regex_type in GENOMIC_DEL_AMBIGUOUS_REGEXPRS:
            match = regex.match(token.change)

            if match:
                match_dict = match.groupdict()
                params.update(match_dict)

                if (
                    classification_type == ClassificationType.GENOMIC_DELETION_AMBIGUOUS
                ):  # noqa: E501
                    self._update_ambiguous_params(params, regex_type)

                    # If ambiguous type not in params, it means we don't support it yet
                    if "ambiguous_type" in params:
                        return GenomicDeletionAmbiguousClassification(**params)
        return None
