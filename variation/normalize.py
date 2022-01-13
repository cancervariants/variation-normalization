"""Module for Variation Normalization."""
from typing import Optional, List
from ga4gh.vrsatile.pydantic.vrsatile_models import VariationDescriptor
from variation.base import Base
from urllib.parse import quote
from variation import logger
from variation.schemas.validation_response_schema import ValidationResult
from variation.schemas.validation_response_schema import ValidationSummary


class Normalize(Base):
    """The Normalize class used to normalize a given variation."""

    @staticmethod
    def get_valid_result(q: str,
                         validations: ValidationSummary,
                         warnings: List) -> ValidationResult:
        """Get valid result from ValidationSummary

        :param str q: Query string
        :param ValidationSummary validations: Validation summary for query
        :param List warnings: List of warnings
        :return: Valid Validation Result
        """
        # For now, only use first valid result
        valid_result = None
        for r in validations.valid_results:
            if r.is_mane_transcript and r.variation:
                valid_result = r
                break
        if not valid_result:
            warning = f"Unable to find MANE Transcript for {q}."
            logger.warning(warning)
            warnings.append(warning)
            valid_result = validations.valid_results[0]
        return valid_result

    def normalize(self, q: str,
                  validations: ValidationSummary,
                  warnings: List) -> Optional[VariationDescriptor]:
        """Normalize a given variation.

        :param str q: The variation to normalize
        :param ValidationSummary validations: Invalid and valid results
        :param List warnings: List of warnings
        :return: An variation descriptor for a valid result if one exists.
            Else, None.
        """
        if not q:
            resp, warnings = self.no_variation_entered()
        else:
            _id = f"normalize.variation:{quote(' '.join(q.strip().split()))}"
            if len(validations.valid_results) > 0:
                valid_result = self.get_valid_result(q, validations, warnings)
                resp, warnings = self.get_variation_descriptor(
                    valid_result.variation, valid_result, _id, warnings)
            else:
                if not q.strip():
                    resp, warnings = self.no_variation_entered()
                else:
                    resp, warnings = self.text_variation_resp(q, _id, warnings)

        self.warnings = warnings
        return resp
