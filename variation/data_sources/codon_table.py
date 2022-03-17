"""Module for Codon Table."""
from typing import Dict, List

from variation.tokenizers.caches import AminoAcidCache


MULTIPLE_CODONS = {"Leu", "Ser", "Arg"}


class CodonTable:
    """Class for codon table data."""

    def __init__(self, amino_acid_cache: AminoAcidCache) -> None:
        """Initialize codon table class.

        :param AminoAcidCache amino_acid_cache: Amino acid code data
        """
        self.table = self._set_codon_table()
        self.amino_acid_cache = amino_acid_cache
        self._dna_to_rna = {
            "T": "A", "A": "U",
            "G": "C", "C": "G"
        }

    @staticmethod
    def _set_codon_table() -> Dict[str, str]:
        """Return RNA codon mappings for amino acids.

        :return: RNA codon --> amino acid single code mappings
        """
        return {
            "AUA": "I", "AUC": "I", "AUU": "I", "AUG": "M",
            "ACA": "T", "ACC": "T", "ACG": "T", "ACU": "T",
            "AAC": "N", "AAU": "N", "AAA": "K", "AAG": "K",
            "AGC": "S", "AGU": "S", "AGA": "R", "AGG": "R",
            "CUA": "L", "CUC": "L", "CUG": "L", "CUU": "L",
            "CCA": "P", "CCC": "P", "CCG": "P", "CCU": "P",
            "CAC": "H", "CAU": "H", "CAA": "Q", "CAG": "Q",
            "CGA": "R", "CGC": "R", "CGG": "R", "CGU": "R",
            "GUA": "V", "GUC": "V", "GUG": "V", "GUU": "V",
            "GCA": "A", "GCC": "A", "GCG": "A", "GCU": "A",
            "GAC": "D", "GAU": "D", "GAA": "E", "GAG": "E",
            "GGA": "G", "GGC": "G", "GGG": "G", "GGU": "G",
            "UCA": "S", "UCC": "S", "UCG": "S", "UCU": "S",
            "UUC": "F", "UUU": "F", "UUA": "L", "UUG": "L",
            "UAC": "Y", "UAU": "Y", "UAA": "*", "UAG": "*",
            "UGC": "C", "UGU": "C", "UGA": "*", "UGG": "W",
        }

    def get_codons(self, amino_acid: str) -> List[str]:
        """Return a list of codons for an amino acid.

        :param str amino_acid: Amino acid to get codons for
        :return: List of codons
        """
        amino_acid = amino_acid.upper()
        if len(amino_acid) == 3:
            amino_acid = self.amino_acid_cache.convert_three_to_one(amino_acid)

        codons = list()
        for k, v in self.table.items():
            if v == amino_acid:
                codons.append(k)
        codons.sort()
        return codons

    def dna_to_rna(self, dna_codon: str) -> str:
        """Convert DNA codon to RNA codon.

        :param str dna_codon: DNA codon
        :return: RNA codon
        """
        dna_codon_list = list(dna_codon)
        rna_codon = ""
        for char in dna_codon_list:
            rna_codon += self._dna_to_rna[char]
        return rna_codon
