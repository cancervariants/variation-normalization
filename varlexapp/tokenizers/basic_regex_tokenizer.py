import re

from .tokenizer import Tokenizer
from ..models import Token

class BasicRegexTokenizer(Tokenizer):
    def __init__(self):
        self.matcher = re.compile(self.pattern(), re.IGNORECASE)

    def match(self, input_string):
        if self.matcher.match(input_string):
            return Token(input_string, self.token_type(), input_string)
        else:
            return None

    def pattern(self):
        pass

    def token_type(self):
        pass
