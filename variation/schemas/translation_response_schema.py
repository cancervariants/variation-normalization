"""Module for Translation Response Schema."""
from typing import Dict, Optional

from cool_seq_tool.schemas import TranscriptPriorityLabel
from pydantic import BaseModel, StrictStr

from variation.schemas.validation_response_schema import ValidationResult


class TranslationResult(BaseModel):
    """Translation Result"""

    vrs_variation: Optional[Dict] = {}
    vrs_seq_loc_ac: Optional[StrictStr] = None
    vrs_seq_loc_ac_status: StrictStr = "na"
    og_ac: Optional[StrictStr] = None
    validation_result: ValidationResult


# Define accession priority. First has highest priority, last has lowest priority
AC_PRIORITY_LABELS = [
    TranscriptPriorityLabel.MANESelect.value,
    TranscriptPriorityLabel.MANEPlusClinical.value,
    TranscriptPriorityLabel.LongestCompatibleRemaining.value,
    "GRCh38",
    "na",
]
