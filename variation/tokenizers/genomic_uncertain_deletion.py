"""A module for tokenizing genomic uncertain deletion."""
from variation.schemas.token_response_schema import \
    GenomicUncertainDeletionToken
from variation.tokenizers.deletion_range_base import DeletionRangeBase


class GenomicUncertainDeletion(DeletionRangeBase):
    """The tokenizer class for genomic uncertain deletion."""

    def _get_parts(self, parts):
        """Set parts for genomic uncertain deletion.

        :param list parts: Parts of input string
        """
        conditions = (
            len(parts) == 4,
            parts[0] == 'g.(?' and parts[3] == '?)del',
            parts[1].endswith(')'),
            parts[2].startswith('(')
        )
        if all(conditions):
            parts[1] = parts[1][:-1]
            parts[2] = parts[2][1:]

            try:
                parts[1] = int(parts[1])
                parts[2] = int(parts[2])
            except ValueError:
                return None
            else:
                if parts[1] < parts[2]:
                    self.parts['start_pos2_del'] = parts[1]
                    self.parts['end_pos1_del'] = parts[2]
        return None

    def return_token(self, params):
        """Return Genomic Uncertain Deletion token."""
        if params['start_pos2_del'] is None or \
                params['end_pos1_del'] is None:
            return None

        if params['reference_sequence'] == 'g':
            params['start_pos1_del'] = "?"
            params['end_pos2_del'] = "?"
            return GenomicUncertainDeletionToken(**params)
