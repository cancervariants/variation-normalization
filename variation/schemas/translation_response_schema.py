"""Module for Translation Response Schema."""
from typing import Optional, Dict

from pydantic import BaseModel, StrictStr
from cool_seq_tool.schemas import TranscriptPriorityLabel


class TranslationResult(BaseModel):
    """Translation Result"""

    vrs_variation: Optional[Dict]
    vrs_seq_loc_ac: Optional[StrictStr]
    vrs_seq_loc_ac_status: StrictStr = "na"


# Used to sort translation results
SORT_AC_ORDER = {
    TranscriptPriorityLabel.MANESelect.value: 0,
    TranscriptPriorityLabel.MANEPlusClinical.value: 1,
    TranscriptPriorityLabel.LongestCompatibleRemaining.value: 2,
    "GRCh38": 3,
    "na": 4
}
