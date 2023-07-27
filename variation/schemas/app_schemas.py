"""Module for schemas used throughout the app"""
from enum import Enum, IntEnum


class Endpoint(str, Enum):
    """Define endpoint names in app that lead to decisions such as hgvs_dup_del_mode
    option.
    """

    TO_VRS = "to_vrs"
    NORMALIZE = "normalize"
    HGVS_TO_COPY_NUMBER_COUNT = "hgvs_to_copy_number_count"
    HGVS_TO_COPY_NUMBER_CHANGE = "hgvs_to_copy_number_change"


class AmbiguousRegexType(IntEnum):
    """Helps determine the regex that was used in ambiguous expressions"""

    REGEX_1 = 1
    REGEX_2 = 2
    REGEX_3 = 3
