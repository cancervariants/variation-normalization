from .basic_regex_tokenizer import BasicRegexTokenizer

class GainOfFunction(BasicRegexTokenizer):
    def pattern(self) -> str:
        return r'\bGAIN[ -]OF[ -]FUNCTION\b'

    def token_type(self) -> str:
        return 'GainOfFunction'
