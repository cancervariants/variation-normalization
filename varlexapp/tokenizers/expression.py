from .basic_regex_tokenizer import BasicRegexTokenizer

class Expression(BasicRegexTokenizer):
    def pattern(self) -> str:
        return r'\b((?<!over )expression|(?<!under )expression)\b'

    def token_type(self) -> str:
        return 'Expression'
