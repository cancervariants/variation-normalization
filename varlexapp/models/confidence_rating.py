from enum import Enum, unique

@unique
class ConfidenceRating(Enum):
    VERY_LOW = 1
    LOW = 2
    MEDIUM = 3
    HIGH = 4

    def __str__(self):
        self.name
