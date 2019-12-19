from enum import Enum, unique

@unique
#not comprehensive yet
class ClassificationType(Enum):
    FUSION = 1
    PROTEIN_SUBSTITUTION = 2
    PROTEIN_FRAMESHIFT = 3
    PROTEIN_ALTERNATE = 4
    PROTEIN_DELINS = 5
    PROTEIN_TERMINATION = 6
    PROTEIN_DUPLICATION = 7
    ONCOGENIC = 8
    EXPRESSION = 9
    COMPLEX = 10

    def __str__(self):
        return self.name
