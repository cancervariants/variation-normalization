"""Module containing schemas used in HGVS To Copy Number endpoints"""
from enum import Enum


class CopyNumberType(str, Enum):
    """Define copy number types"""

    RELATIVE = "relative_copy_number"
    ABSOLUTE = "absolute_copy_number"
