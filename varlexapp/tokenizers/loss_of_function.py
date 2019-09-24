from .basic_regex_tokenizer import BasicRegexTokenizer

class LossOfFunction(BasicRegexTokenizer):
    def pattern(self):
        return r'\bLOSS[ -]OF[ -]FUNCTION\b'

    def token_type(self):
        return 'LossOfFunction'
