from typing import Set

class AminoAcidCache:
    def __init__(self) -> None:
        self.__amino_acid_codes = self.__load_amino_acid_codes()

    def __contains__(self, item: str) -> bool:
        return item in self.__amino_acid_codes

    # TODO: Use BioPython
    def __load_amino_acid_codes(self) -> Set[str]:
        return ({'Ala', 'Asx', 'Cys', 'Asp', 'Glu', 'Phe', 'Gly', 'His', 'Ile',
                 'Lys', 'Leu', 'Met', 'Asn', 'Pro', 'Gln', 'Arg', 'Ser', 'Thr',
                 'Sec', 'Val', 'Trp', 'Xaa', 'Tyr', 'Glx'} -
                {'Asx', 'Sec', 'Glx', 'Xaa'}) | set('ACDEFGHIKLMNPQRSTVWY')
