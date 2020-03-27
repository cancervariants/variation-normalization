from typing import Set

from Bio.Alphabet import IUPAC, ThreeLetterProtein

class AminoAcidCache:
    def __init__(self) -> None:
        self.__amino_acid_codes = self.__load_amino_acid_codes()

    def __contains__(self, item: str) -> bool:
        return item in self.__amino_acid_codes

    def __load_amino_acid_codes(self) -> Set[str]:
        return (set(ThreeLetterProtein.letters) - { 'Asx', 'Sec', 'Glx', 'Xaa' }) | set(IUPAC.IUPACProtein.letters)
