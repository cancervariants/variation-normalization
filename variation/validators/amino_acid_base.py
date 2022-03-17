"""Module for commonly used validator methods for protein references."""
from typing import List

from variation.tokenizers.caches import AminoAcidCache
from variation.data_sources import SeqRepoAccess


class AminoAcidBase:
    """Amino Acid Base class."""

    def __init__(self, seq_repo_access: SeqRepoAccess,
                 amino_acid_cache: AminoAcidCache) -> None:
        """Initialize Amino Acid Base class.

        :param SeqRepoAccess seq_repo_access: Access to SeqRepo data
        :param AminoAcidCache amino_acid_cache: Amino acid code data
        """
        self.seqrepo_access = seq_repo_access
        self.amino_acid_cache = amino_acid_cache

    def check_ref_aa(self, t: str, aa: str, pos: int, errors: List) -> None:
        """Check that reference amino acid matches actual amino acid.

        :param str t: Transcript
        :param str aa: Expected Amino Acid
        :param int pos: Expected position
        :param List errors: List of errors
        """
        ref_aa = self.seqrepo_access.get_sequence(t, pos)
        if ref_aa and len(ref_aa) == 1 \
                and len(aa) == 3:
            aa = self.amino_acid_cache.convert_three_to_one(aa)

        if ref_aa != aa:
            errors.append(f"Needed to find {aa} at "
                          f"position {pos} on {t} "
                          f"but found {ref_aa}")
