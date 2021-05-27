"""Module for Codon Table."""
from typing import Dict, List

from variant.tokenizers.caches import AminoAcidCache


class CodonTable:
    """Class for codon table data."""

    def __init__(self):
        """Initialize codon table class."""
        self.codon_table = self._set_codon_table()
        self.amino_acid_cache = AminoAcidCache()

    @staticmethod
    def _set_codon_table() -> Dict[str, str]:
        """Return codon mappings for amino acids.

        :return: codon --> amino acid single code mappings
        """
        return {
            'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
            'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
        }

    def _get_codons(self, amino_acid) -> List[str]:
        """Return a list of codons for an amino acid.

        :param str amino_acid: Amino acid to get codons for
        :return: List of codons
        """
        amino_acid = amino_acid.upper()
        if len(amino_acid) == 3:
            amino_acid = self.amino_acid_cache.convert_three_to_one(amino_acid)

        codons = list()
        for k, v in self.codon_table.items():
            if v == amino_acid:
                codons.append(k)
        codons.sort()
        return codons
