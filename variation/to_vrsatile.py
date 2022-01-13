"""Module for toVRSATILE endpoint"""
from typing import List, Tuple
from urllib.parse import quote
from ga4gh.vrsatile.pydantic.vrsatile_models import Extension
from variation.schemas.validation_response_schema import ValidationResult
from variation.to_vrs import ToVRS
from variation import logger
from variation.schemas.validation_response_schema import ValidationSummary
from variation.schemas.normalize_response_schema import ReferenceGenomeBuild


class ToVRSATILE(ToVRS):
    """ToVRSATILE class"""

    @staticmethod
    def get_valid_result(
            q: str, validations: ValidationSummary,
            warnings: List) -> Tuple[ValidationResult, List, List]:
        """Get valid results
        :param str q: Query string
        :param ValidationSummary validations: Validation Summary for query
        :param List warnings: List of warnings
        :return: Valid result, list of valid results, list of accession status
        """
        list_of_valid_res = list()

        for r in validations.valid_results:
            if r.is_valid and r.variation:
                list_of_valid_res.append(r)

        if not list_of_valid_res:
            warning = f"Unable to normalize variants for {q}."
            logger.warning(warning)
            warnings.append(warning)
            list_of_valid_res = [validations.valid_results[0]]
        return list_of_valid_res

    def to_vrsatile(self, q: str, validations: ValidationSummary,
                    warnings: List) -> Tuple[List, List]:
        """Get valid results

        :param str q: Query string
        :param ValidationSummary validations: Validation Summary for query
        :param List warnings: List of warnings
        :return: Response, list of warnings
        """
        resp = list()
        if not q:
            resp, warnings = self.no_variation_entered()
            resp = [resp]
        else:
            _id = f"normalize.variation:{quote(' '.join(q.strip().split()))}"
            if len(validations.valid_results) > 0:
                list_of_valid_res = \
                    self.get_valid_result(q, validations, warnings)

                for val_res in list_of_valid_res:
                    tovrsatile_resp, _ = \
                        self.get_variation_descriptor(
                            val_res.variation, val_res, _id,
                            warnings)

                    # we need a better answer for reference genome build
                    # because all the possible accessions are from transcript mapping file  # noqa: E501
                    # if users chose to use GRCh37, all the possible accessions will be in GRCh37  # noqa: E501
                    # but the accession_id for mane transcript is a GRCh38
                    # it would be confusing
                    # also need to consider adding mane_transcript ENST_ id to variation descriptor  # noqa: E501
                    if val_res.is_mane_transcript:
                        ref_genome = ReferenceGenomeBuild.GRCH38
                        mane_status = "mane transcript"
                    else:
                        ref_genome = ReferenceGenomeBuild.GRCH37
                        mane_status = "N/A"
                    tovrsatile_resp.extensions = [
                        Extension(
                            name='possible accessions',
                            value=val_res.possible_ac
                        ).dict(exclude_none=True),
                        Extension(
                            name='human reference genome assembly',
                            value=ref_genome
                        ).dict(exclude_none=True),
                        Extension(
                            name='mane status',
                            value=mane_status
                        ).dict(exclude_none=True)
                    ]
                    resp.append(tovrsatile_resp)
            else:
                if not q.strip():
                    resp, warnings = self.no_variation_entered()
                    resp = [resp]
                else:
                    resp, _ = self.text_variation_resp(
                        q, _id, warnings)
                    resp = [resp]

        self.warnings = warnings
        return resp, warnings
