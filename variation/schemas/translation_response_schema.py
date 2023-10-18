"""Module for Translation Response Schema."""
from enum import Enum
from typing import Dict, Optional

from cool_seq_tool.schemas import TranscriptPriority
from pydantic import BaseModel, StrictStr

from variation.schemas.validation_response_schema import ValidationResult


class VrsSeqLocAcStatus(str, Enum):
    """Create enum for VRS SequenceLocation accession status.
    Order when defining matters.
        First has highest priority, last has lowest priority
    """

    MANE_SELECT = TranscriptPriority.MANE_SELECT.value
    MANE_PLUS_CLINICAL = TranscriptPriority.MANE_PLUS_CLINICAL.value
    LONGEST_COMPATIBLE_REMAINING = TranscriptPriority.LONGEST_COMPATIBLE_REMAINING.value
    GRCH38 = TranscriptPriority.GRCH38.value
    NA = "na"


AC_PRIORITY_LABELS = [m for m in VrsSeqLocAcStatus.__members__.values()]


class TranslationResult(BaseModel):
    """Translation Result"""

    vrs_variation: Optional[Dict] = None
    vrs_seq_loc_ac: Optional[StrictStr] = None
    vrs_seq_loc_ac_status: VrsSeqLocAcStatus = VrsSeqLocAcStatus.NA
    og_ac: Optional[StrictStr] = None
    validation_result: ValidationResult
