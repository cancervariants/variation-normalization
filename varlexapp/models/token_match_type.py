from enum import Enum, unique

@unique
class TokenMatchType(Enum):
    ID = 1
    SYMBOL = 2
    ALIAS = 3
    PREVIOUS = 4
    UNSPECIFIED = 5


    def __str__(self):
        self.name
