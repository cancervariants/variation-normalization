"""Module for Translation Response Schema."""
from typing import Optional, Dict

from pydantic import BaseModel


class TranslationResult(BaseModel):
    """Translation Result"""

    vrs_variation: Optional[Dict]
    vrs_seq_loc_ac: Optional[str]
