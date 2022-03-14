"""Module containing schemas used in HGVS To Copy Number endpoints"""
from enum import Enum


class CopyNumberType(str, Enum):
    """Define copy number types"""

    RELATIVE = "relative_copy_number"
    ABSOLUTE = "absolute_copy_number"


class RelativeCopyClass(str, Enum):
    """The relative copy class"""

    COMPLETE_LOSS = "complete loss"
    PARTIAL_LOSS = "partial loss"
    COPY_NEUTRAL = "copy neutral"
    LOW_LEVEL_GAIN = "low-level gain"
    HIGH_LEVEL_GAIN = "high-level gain"
