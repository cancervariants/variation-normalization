"""Module for Translation Response Schema."""
from enum import Enum
from typing import Dict, Optional

from cool_seq_tool.schemas import TranscriptPriorityLabel
from pydantic import BaseModel, StrictStr

from variation.schemas.validation_response_schema import ValidationResult


class VrsSeqLocAcStatus(str, Enum):
    """Create enum for VRS SequenceLocation accession status.
    Order when defining matters.
        First has highest priority, last has lowest priority
    Once issue-191 is resolved in cool-seq-tool, we should use the
    TranscriptPriorityLabel enum
    """

    MANE_SELECT = TranscriptPriorityLabel.MANESelect.value
    MANE_PLUS_CLINICAL = TranscriptPriorityLabel.MANEPlusClinical.value
    LONGEST_COMPATIBLE_REMAINING = (
        TranscriptPriorityLabel.LongestCompatibleRemaining.value
    )
    GRCH38 = "GRCh38"  # will change to lowercase in cool-seq-tool issue-191
    NA = "na"


AC_PRIORITY_LABELS = [m for m in VrsSeqLocAcStatus.__members__.values()]


class TranslationResult(BaseModel):
    """Translation Result"""

    vrs_variation: Optional[Dict]
    vrs_seq_loc_ac: Optional[StrictStr]
    vrs_seq_loc_ac_status: VrsSeqLocAcStatus = VrsSeqLocAcStatus.NA
    og_ac: Optional[StrictStr]
    validation_result: ValidationResult
