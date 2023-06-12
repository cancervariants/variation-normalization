"""A module for the HGVS Classifier."""
from typing import List, Optional

from variation.schemas.classification_response_schema import (
    Classification, ClassificationType, GenomicSubstitutionClassification,
    GenomicInsertionClassification, ProteinSubstitutionClassification,
    CdnaSubstitutionClassification, CdnaDeletionClassification,
    ProteinDeletionClassification, ProteinDelInsClassification,
    ProteinInsertionClassification, ProteinReferenceAgreeClassification,
    CdnaDelInsClassification, GenomicDelInsClassification, CdnaInsertionClassification,
    Nomenclature, SequenceOntology
)
from variation.schemas.token_response_schema import HgvsToken, TokenType, CoordinateType
from variation.regex import GENOMIC_REGEXPRS, PROTEIN_REGEXPRS, CDNA_REGEXPRS
from variation.classifiers import Classifier


class HgvsClassifier(Classifier):
    """The Genomic Substitution Classifier class."""

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
        elif token.coordinate_type == CoordinateType.CODING_DNA:
            classification = self._cdna_classification(token, params)
        elif token.coordinate_type == CoordinateType.PROTEIN:
            classification = self._protein_classification(token, params)

        return classification

    def _protein_classification(self, token, params) -> Optional[Classification]:
        for regex, _, classification_type in PROTEIN_REGEXPRS:
            match = regex.match(token.change)

            if match:
                match_dict = match.groupdict()
                params.update(match_dict)

                if classification_type == ClassificationType.PROTEIN_SUBSTITUTION:
                    params["pos"] = int(params["pos"])
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

    def _cdna_classification(self, token, params) -> Optional[Classification]:
        for regex, _, classification_type in CDNA_REGEXPRS:
            match = regex.match(token.change)

            if match:
                match_dict = match.groupdict()
                params.update(match_dict)

                if classification_type == ClassificationType.CODING_DNA_SUBSTITUTION:
                    len_ref = params["ref"]
                    len_alt = params["alt"]

                    if len_ref == 1 and len_alt == 1:
                        params["so_id"] = SequenceOntology.SNV
                    else:
                        params["so_id"] = SequenceOntology.MNV

                    params["pos"] = int(params["pos"])

                    return CdnaSubstitutionClassification(**params)
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

    def _genomic_classification(self, token, params) -> Optional[Classification]:
        for regex, _, classification_type in GENOMIC_REGEXPRS:
            match = regex.match(token.change)

            if match:
                match_dict = match.groupdict()
                params.update(match_dict)

                if classification_type == ClassificationType.GENOMIC_SUBSTITUTION:
                    len_ref = params["ref"]
                    len_alt = params["alt"]

                    if len_ref == 1 and len_alt == 1:
                        params["so_id"] = SequenceOntology.SNV
                    else:
                        params["so_id"] = SequenceOntology.MNV

                    params["pos"] = int(params["pos"])

                    return GenomicSubstitutionClassification(**params)
                elif classification_type == ClassificationType.GENOMIC_DELINS:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = int(params["pos1"]) if params["pos1"] is not None else params["pos1"]  # noqa: E501
                    return GenomicDelInsClassification(**params)
                elif classification_type == ClassificationType.GENOMIC_INSERTION:
                    params["pos0"] = int(params["pos0"])
                    params["pos1"] = int(params["pos1"]) if params["pos1"] is not None else params["pos1"]  # noqa: E501
                    return GenomicInsertionClassification(**params)

        return None
