"""Module for schemas used throughout the app"""
from enum import Enum


class Endpoint(str, Enum):
    """Define contrains for endpoint names"""

    TO_VRS = "to_vrs"
    NORMALIZE = "normalize"
    TRANSLATE_IDENTIFIER = "translate_identifier"
    GNOMAD_VCF_TO_PROTEIN = "gnomad_vcf_to_protein"
    TRANSLATE_FROM = "translate_from"
    HGVS_TO_COPY_NUMBER_COUNT = "hgvs_to_copy_number_count"
    HGVS_TO_COPY_NUMBER_CHANGE = "hgvs_to_copy_number_change"
    TO_CANONICAL = "to_canonical"
