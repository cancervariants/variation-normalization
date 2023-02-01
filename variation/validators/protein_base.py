"""Module for commonly used validator methods for protein references."""
from typing import List

from cool_seq_tool.data_sources import SeqRepoAccess
from bioutils.sequences import aa3_to_aa1


class ProteinBase:
    """Protein Base class."""

    def __init__(self, seq_repo_access: SeqRepoAccess) -> None:
        """Initialize Protein Base class.

        :param SeqRepoAccess seq_repo_access: Access to SeqRepo data
        """
        self.seqrepo_access = seq_repo_access

    def check_ref_aa(self, t: str, aa: str, pos: int, errors: List) -> None:
        """Check that reference protein matches actual protein.

        :param str t: Transcript
        :param str aa: Expected Protein
        :param int pos: Expected position
        :param List errors: List of errors
        """
        ref_aa, warning = self.seqrepo_access.get_reference_sequence(t, pos)
        if not ref_aa:
            errors.append(warning)
        else:
            if ref_aa and len(ref_aa) == 1 and len(aa) == 3:
                aa = aa3_to_aa1(aa.capitalize())

            if ref_aa != aa:
                errors.append(f"Needed to find {aa} at position {pos} on {t} but found"
                              f" {ref_aa}")
