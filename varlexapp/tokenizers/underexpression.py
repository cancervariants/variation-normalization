from .basic_regex_tokenizer import BasicRegexTokenizer

class UnderExpression(BasicRegexTokenizer):
    def pattern(self) -> str:
        return r'\bunder ?expression\b'

    def token_type(self) -> str:
        return 'UnderExpression'
