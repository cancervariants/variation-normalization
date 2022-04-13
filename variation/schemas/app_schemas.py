"""Module for schemas used throughout the app"""
from enum import Enum


class Endpoint(str, Enum):
    """Define contrains for endpoint names"""

    TO_VRS = "to_vrs"
    NORMALIZE = "normalize"
    TRANSLATE_IDENTIFIER = "translate_identifier"
    GNOMAD_VCF_TO_PROTEIN = "gnomad_vcf_to_protein"
    TRANSLATE_FROM = "translate_from"
    HGVS_TO_ABSOLUTE_CN = "hgvs_to_absolute_copy_number"
    HGVS_TO_RELATIVE_CN = "hgvs_to_relative_copy_number"
    TO_CANONICAL = "to_canonical"
