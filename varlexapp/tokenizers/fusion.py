from .basic_regex_tokenizer import BasicRegexTokenizer

class Fusion(BasicRegexTokenizer):
    def pattern(self) -> str:
        return r'\bfusion(s)?\b'

    def token_type(self) -> str:
        return 'Fusion'
