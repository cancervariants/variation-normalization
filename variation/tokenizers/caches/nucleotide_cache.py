"""A module to cache nucleotides."""


class NucleotideCache:
    """A class to cache nucleotides."""

    def __init__(self) -> None:
        """Initialize NucleotideCache class.
        https://varnomen.hgvs.org/bg-material/standards/
        """
        self.base_nucleotides = {"A", "C", "T", "G"}
        self.nucleotides = {
            "A": ["A"],
            "C": ["C"],
            "T": ["T"],
            "G": ["G"],
            "B": ["C", "G", "T"],
            "D": ["A", "G", "T"],
            "H": ["A", "C", "T"],
            "K": ["G", "T"],
            "M": ["A", "C"],
            "N": ["A", "C", "G", "T"],
            "R": ["A", "G"],
            "S": ["G", "C"],
            "V": ["A", "C", "G"],
            "W": ["A", "T"],
            "Y": ["C", "T"],
            "X*": ["A", "C", "G", "T"],
            "-*": []
        }
