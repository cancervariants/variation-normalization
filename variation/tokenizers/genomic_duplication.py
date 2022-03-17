"""A module for Genomic Duplication Tokenization."""
from typing import Optional, Union, List

from variation.schemas.token_response_schema import DuplicationAltType, \
    GenomicDuplicationToken, GenomicDuplicationRangeToken
from variation.tokenizers.duplication_base import DuplicationBase


class GenomicDuplication(DuplicationBase):
    """Class for tokenizing duplications on the genomic coordinate."""

    def _get_parts(self, parts: List) -> None:
        if len(parts) != 1 or not parts[0].startswith("g."):
            return None

        parts[0] = parts[0][2:]
        if "_" in parts[0]:
            if parts[0].count("_") == 1:
                pos = parts[0].split("_")
                try:
                    pos[0] = int(pos[0])
                    pos[1] = int(pos[1])
                except ValueError:
                    pass
                else:
                    if pos[0] < pos[1]:
                        self.parts["start_pos1_dup"] = pos[0]
                        self.parts["start_pos2_dup"] = pos[1]
                        self.parts["coordinate_type"] = "g"
            else:
                self.parts["alt_type"] = DuplicationAltType.DUPLICATION_RANGE
                parts = parts[0].split("_")
                len_parts = len(parts)
                if len_parts == 4:
                    for part_ix, parts_field in [
                        (0, "start_pos1_dup"),
                        (1, "start_pos2_dup"),
                        (2, "end_pos1_dup"),
                        (3, "end_pos2_dup")
                    ]:
                        part_val = self._check_uncertain_or_int(parts[part_ix])
                        if part_val is None:
                            return None
                        else:
                            self.parts[parts_field] = part_val
                    self.parts["coordinate_type"] = "g"

                elif len_parts == 3:
                    if "(" in parts[0] and ")" in parts[1]:
                        # Format is: (?_#)_#dup
                        for part_ix, parts_field in [
                            (0, "start_pos1_dup"),
                            (1, "start_pos2_dup"),
                            (2, "end_pos1_dup")
                        ]:
                            part_val = self._check_uncertain_or_int(
                                parts[part_ix])
                            if part_val is None:
                                return None
                            else:
                                self.parts[parts_field] = part_val
                    else:
                        # Format is #_(#_?)dup
                        for part_ix, parts_field in [
                            (0, "start_pos1_dup"),
                            (1, "end_pos1_dup"),
                            (2, "end_pos2_dup")
                        ]:
                            part_val = self._check_uncertain_or_int(
                                parts[part_ix])
                            if part_val is None:
                                return None
                            else:
                                self.parts[parts_field] = part_val
                    self.parts["coordinate_type"] = "g"
        else:
            try:
                pos = int(parts[0])
            except ValueError:
                pass
            else:
                self.parts["start_pos1_dup"] = pos
                self.parts["coordinate_type"] = "g"
        return None

    def _check_uncertain_or_int(self, part: str) -> Optional[Union[int, str]]:
        part = part.replace("(", "")
        part = part.replace(")", "")
        try:
            return int(part)
        except ValueError:
            if part == "?":
                self.parts["alt_type"] = \
                    DuplicationAltType.UNCERTAIN_DUPLICATION
                return part
        return None

    def return_token(self) -> Optional[Union[GenomicDuplicationRangeToken,
                                             GenomicDuplicationToken]]:
        """Return token instance if a match is found."""
        # we only set this field if it"s valid
        if self.parts["coordinate_type"] == "g":
            if self.parts["end_pos1_dup"] is None and \
                    self.parts["end_pos2_dup"] is None:

                if self.parts["start_pos2_dup"] and self.parts["start_pos1_dup"] > self.parts["start_pos2_dup"]:  # noqa: E501
                    return

                self.parts["alt_type"] = DuplicationAltType.DUPLICATION
                return GenomicDuplicationToken(**self.parts)
            else:
                prev_val = None
                for field in ["start_pos1_dup", "start_pos2_dup",
                              "end_pos1_dup", "end_pos2_dup"]:
                    val = self.parts[field]
                    if val not in ["?", None]:
                        if prev_val is not None:
                            if prev_val > val:
                                return
                        prev_val = val
                return GenomicDuplicationRangeToken(**self.parts)
