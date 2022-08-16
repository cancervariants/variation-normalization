"""Module for general functionality throughout the app"""
from typing import List, Tuple, Optional

from ga4gh.vrsatile.pydantic.vrs_models import Text
from ga4gh.vrsatile.pydantic.vrsatile_models import VariationDescriptor
from ga4gh.core import ga4gh_identify
from ga4gh.vrs import models

from variation.schemas.validation_response_schema import ValidationSummary, \
    ValidationResult
from variation.schemas.token_response_schema import Token


def no_variation_entered() -> Tuple[None, List[str]]:
    """Return response when no variation queried.

    :return: None, list of warnings
    """
    warnings = ["No variation was entered to normalize"]
    return None, warnings


def get_mane_valid_result(q: str, validations: ValidationSummary,
                          warnings: List) -> ValidationResult:
    """Get valid result from ValidationSummary

    :param str q: Query string
    :param ValidationSummary validations: Validation summary for query
    :param List warnings: List of warnings
    :return: Valid Validation Result
    """
    # For now, only use first valid result
    valid_result = None
    if validations and validations.valid_results:
        for r in validations.valid_results:
            if r.is_mane_transcript and r.variation:
                valid_result = r
                break
    if not valid_result:
        warning = f"Unable to find MANE Transcript for {q}."
        warnings.append(warning)
        if validations and validations.valid_results:
            valid_result = validations.valid_results[0]
    return valid_result


def no_variation_resp(
    label: str, _id: str, warnings: List, untranslatable_returns_text: bool = False
) -> Tuple[Optional[VariationDescriptor], List]:
    """Return Variation Descriptor with variation set as Text or return `None` for
    queries that could not be normalized

    :param str label: Initial input query
    :param str _id: _id field for variation descriptor
    :param List warnings: List of warnings
    :param bool untranslatable_returns_text: `True` return VRS Text Object when
        unable to translate or normalize query. `False` return `None` when
        unable to translate or normalize query.
    :return: Variation descriptor or `None`, warnings
    """
    warning = f"Unable to translate {label}"
    if untranslatable_returns_text:
        text = models.Text(definition=label, type="Text")
        text._id = ga4gh_identify(text)
        variation = Text(**text.as_dict())
        resp = VariationDescriptor(id=_id, label=label, variation=variation)
    else:
        resp = None

    if not warnings:
        warnings.append(warning)
    return resp, warnings


def is_token_type(valid_result_tokens: List, token_type: str) -> bool:
    """Return whether or not token_type is in valid_result_tokens.

    :param List valid_result_tokens: Valid token matches
    :param str token_type: The token's type
    :return: Whether or not token_type is in valid_result_tokens
    """
    for t in valid_result_tokens:
        if t.token_type == token_type:
            return True
    return False


def get_instance_type_token(valid_result_tokens: List,
                            instance_type: Token) -> Optional[Token]:
    """Return the tokens for a given instance type.

    :param List valid_result_tokens: A list of valid tokens for the input
        string
    :param Token instance_type: The instance type to check
    :return: Token for a given instance type
    """
    for t in valid_result_tokens:
        if isinstance(t, instance_type):
            return t
    return None
