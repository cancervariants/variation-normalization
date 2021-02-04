from .token import Token


class ProteinSubstitutionToken(Token):
    def __init__(self, input_string: str, ref_protein: str, alt_protein: str,
                 pos: int) -> None:
        self.ref_protein = ref_protein
        self.alt_protein = alt_protein
        self.pos = pos
        super().__init__(input_string, 'ProteinSubstitution', input_string)
