"""A module for creating a Protein Substitution Token instance."""
from .token import Token


class ProteinSubstitutionToken(Token):
    """The ProteinSubstitutionToken class."""

    def __init__(self, input_string: str, ref_protein: str, alt_protein: str,
                 pos: int) -> None:
        """Initialize the ProteinSubstitutionToken class.

        :param str input_string: The input string
        :param str ref_protein: The reference amino acid
        :param str alt_protein: The new amino acid
        :param int pos: The position of the amino acid substituted
        """
        self.ref_protein = ref_protein
        self.alt_protein = alt_protein
        self.pos = pos
        super().__init__(input_string, 'ProteinSubstitution', input_string)
