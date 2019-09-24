from Bio.Alphabet import IUPAC, ThreeLetterProtein

class AminoAcidCache:
    def __init__(self):
        self.__amino_acid_codes = self.__load_amino_acid_codes()

    def __contains__(self, item):
        return item in self.__amino_acid_codes

    def __load_amino_acid_codes(self):
        return (set(ThreeLetterProtein.letters) - { 'Asx', 'Sec', 'Glx', 'Xaa' }) | set(IUPAC.IUPACProtein.letters)
