from enum import Enum, unique

@unique
class ConfidenceRating(Enum):
    INTERSECTION = 1
    SUPERSET = 2
    UNORDERED = 3
    EXACT = 4

    def __str__(self):
        self.name
