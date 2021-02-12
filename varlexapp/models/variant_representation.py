"""Module for Variant Representation model."""
from .location import Location
from .sequence_state import SequenceState


class VariantRepresentation:
    """The Variant Representation model class."""

    def __init__(self, _id: str, location: Location, state: SequenceState,
                 vrs_type: str) -> None:
        """Initialize variant representation."""
        self._id = _id
        self.state = state
        self.location = location
        self.type = vrs_type
