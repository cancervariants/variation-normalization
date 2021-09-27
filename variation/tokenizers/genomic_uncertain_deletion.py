"""A module for tokenizing genomic uncertain deletion."""
from pydantic.error_wrappers import ValidationError

from variation.schemas.token_response_schema import \
    GenomicUncertainDeletionToken
from variation.tokenizers.deletion_range_base import DeletionRangeBase


class GenomicUncertainDeletion(DeletionRangeBase):
    """The tokenizer class for genomic uncertain deletion."""

    def _get_parts(self, parts):
        """Set parts for genomic uncertain deletion.

        :param list parts: Parts of input string
        """
        len_parts = len(parts)
        if len_parts not in [3, 4]:
            return None

        if not parts[0].startswith('g.'):
            return None

        parts[0] = parts[0][2:]
        parts[len_parts - 1] = parts[len_parts - 1][:-3]

        if len_parts == 3:
            if '(' in parts[0] and ')' in parts[1]:
                # Format is: (?_#)_#del
                for part_ix, parts_field in [
                    (0, 'start_pos1_del'),
                    (1, 'start_pos2_del'),
                    (2, 'end_pos1_del')
                ]:
                    part_val = self._check_uncertain_or_int(
                        parts[part_ix])
                    if part_val is None:
                        return None
                    else:
                        self.parts[parts_field] = part_val
            else:
                # Format is #_(#_?)del
                for part_ix, parts_field in [
                    (0, 'start_pos1_del'),
                    (1, 'end_pos1_del'),
                    (2, 'end_pos2_del')
                ]:
                    part_val = self._check_uncertain_or_int(
                        parts[part_ix])
                    if part_val is None:
                        return None
                    else:
                        self.parts[parts_field] = part_val
            self.parts['reference_sequence'] = 'g'
        elif len_parts == 4:
            for part_ix, parts_field in [
                (0, 'start_pos1_del'),
                (1, 'start_pos2_del'),
                (2, 'end_pos1_del'),
                (3, 'end_pos2_del')
            ]:
                part_val = self._check_uncertain_or_int(parts[part_ix])
                if part_val is None:
                    return None
                else:
                    self.parts[parts_field] = part_val
            self.parts['reference_sequence'] = 'g'
        return None

    def _check_uncertain_or_int(self, part):
        part = part.replace('(', '')
        part = part.replace(')', '')
        try:
            return int(part)
        except ValueError:
            if part == '?':
                return part
        return None

    def return_token(self, params):
        """Return Genomic Uncertain Deletion token."""
        if params['reference_sequence'] == 'g':
            conditions = (
                isinstance(params['start_pos1_del'], int),
                isinstance(params['start_pos2_del'], int),
                isinstance(params['end_pos1_del'], int),
                isinstance(params['end_pos1_del'], int)
            )
            if not all(conditions):
                try:
                    return GenomicUncertainDeletionToken(**params)
                except ValidationError:
                    return None
