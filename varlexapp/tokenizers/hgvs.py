from .tokenizer import Tokenizer
from ..models import Token

from hgvs.parser import Parser
from hgvs.exceptions import HGVSParseError, HGVSInvalidVariantError
from hgvs.validator import IntrinsicValidator

class HGVS(Tokenizer):
    def __init__(self):
        self.parser = Parser()
        self.validator = IntrinsicValidator()

    def match(self, input_string):
        try:
            variant = self.parser.parse_hgvs_variant(input_string)
            self.validator.validate(variant)
            return Token(input_string, 'HGVS')
        except (HGVSParseError, HGVSInvalidVariantError):
            return None
