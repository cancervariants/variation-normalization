"""Module for commonly used validator methods for protein references."""
from variant.tokenizers.caches import AminoAcidCache
from variant.data_sources import SeqRepoAccess


class AminoAcidBase:
    """Amino Acid Base class."""

    def __init__(self, seq_repo_access: SeqRepoAccess,
                 amino_acid_cache: AminoAcidCache):
        """Initialize Amino Acid Base class.

        :param SeqRepoAccess seq_repo_access: Access to SeqRepo data
        :param AminoAcidCache amino_acid_cache: Amino acid code data
        """
        self.seqrepo_access = seq_repo_access
        self.amino_acid_cache = amino_acid_cache

    def check_ref_aa(self, t, aa, pos, errors):
        """Check that reference amino acid matches actual amino acid.

        :param string t: Transcript
        :param str aa: Expected Amino Acid
        :param str pos: Expected position
        :param list errors: List of errors
        """
        ref_aa = self.seqrepo_access.sequence_at_position(t, pos)
        if ref_aa and len(ref_aa) == 1 \
                and len(aa) == 3:
            ref_aa = self.amino_acid_cache.amino_acid_code_conversion[ref_aa]

        if ref_aa != aa:
            errors.append(f"Needed to find {aa} at "
                          f"position {pos} on {t} "
                          f"but found {ref_aa}")
