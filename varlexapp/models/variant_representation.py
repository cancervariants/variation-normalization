from . import ValidationResult, ClassificationType
from .location import Location
from .sequence_state import SequenceState
from ..data_sources import SeqRepoAccess


class VariantRepresentation:
    def __init__(self, location: Location, state: SequenceState) -> None:
        self.state = state
        self.location = location
        self.type = 'Allele'

