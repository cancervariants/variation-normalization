from .basic_regex_tokenizer import BasicRegexTokenizer

class LossOfFunction(BasicRegexTokenizer):
    def pattern(self) -> str:
        return r'\bLOSS[ -]OF[ -]FUNCTION\b'

    def token_type(self) -> str:
        return 'LossOfFunction'
