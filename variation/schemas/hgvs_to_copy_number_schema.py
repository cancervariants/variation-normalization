"""Module containing schemas used in HGVS To Copy Number endpoints"""
from enum import Enum

from ga4gh.vrsatile.pydantic.vrs_models import RelativeCopyClass

from variation.schemas.classification_response_schema import ClassificationType


VALID_CLASSIFICATION_TYPES = [
    ClassificationType.GENOMIC_DUPLICATION,
    ClassificationType.GENOMIC_DELETION,
    ClassificationType.GENOMIC_DELETION_RANGE,
    ClassificationType.GENOMIC_UNCERTAIN_DELETION
]


class CopyNumberType(str, Enum):
    """Define copy number types"""

    RELATIVE = "relative_copy_number"
    ABSOLUTE = "absolute_copy_number"


VALID_RELATIVE_COPY_CLASS = [rcc.value for
                             rcc in RelativeCopyClass.__members__.values()]
