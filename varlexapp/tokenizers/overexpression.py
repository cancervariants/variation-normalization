from .basic_regex_tokenizer import BasicRegexTokenizer

class OverExpression(BasicRegexTokenizer):
    def pattern(self) -> str:
        return r'\bover ?expression\b'

    def token_type(self) -> str:
        return 'OverExpression'
