from .simple_interval import SimpleInterval

class Location:
    def __init__(self, sequence_id: str, interval: SimpleInterval) -> None:
        self.sequence_id = sequence_id
        self.interval = interval
        self.type = 'SequenceLocation'
