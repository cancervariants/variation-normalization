from enum import Enum, unique

@unique
class LookupType(Enum):
    GENE_SYMBOL = 1

    def __str__(self) -> str:
        return self.name
