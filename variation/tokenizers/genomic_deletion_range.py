"""A module for Genomic Deletion Range Tokenization."""
from typing import Dict, Optional, List

from pydantic import NonNegativeFloat

from variation.schemas.token_response_schema import GenomicDeletionRangeToken
from variation.tokenizers.deletion_range_base import DeletionRangeBase


class GenomicDeletionRange(DeletionRangeBase):
    """Class for tokenizing deletion range at the genomic coordinate."""

    def _get_parts(self, parts: List) -> NonNegativeFloat:
        """Set parts for genomic deletion range.

        :param List parts: Parts of input string
        """
        if len(parts) != 4:
            return None

        conditions = (
            parts[0].startswith("g.("),
            parts[1].endswith(")"),
            parts[2].startswith("("),
            parts[3].endswith(")del")
        )

        if all(conditions):
            parts[0] = parts[0][3:]
            parts[1] = parts[1][:-1]
            parts[2] = parts[2][1:]
            parts[3] = parts[3][:-4]

            try:
                parts[0] = int(parts[0])
                parts[1] = int(parts[1])
                parts[2] = int(parts[2])
                parts[3] = int(parts[3])
            except ValueError:
                return None
            else:
                prev_val = None
                for i in range(4):
                    val = parts[i]
                    if val not in ["?", None]:
                        if prev_val is not None:
                            if prev_val > val:
                                return None
                    prev_val = val

                self.parts["start_pos1_del"] = parts[0]
                self.parts["start_pos2_del"] = parts[1]
                self.parts["end_pos1_del"] = parts[2]
                self.parts["end_pos2_del"] = parts[3]
        return None

    def return_token(self, params: Dict) -> Optional[GenomicDeletionRangeToken]:
        """Return Genomic Deletion Range token."""
        conditions = (
            params["start_pos1_del"] is not None,
            params["start_pos2_del"] is not None,
            params["end_pos1_del"] is not None,
            params["end_pos2_del"] is not None
        )
        if all(conditions):
            if params["coordinate_type"] == "g":
                return GenomicDeletionRangeToken(**params)
