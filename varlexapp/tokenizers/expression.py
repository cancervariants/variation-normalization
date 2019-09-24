from .basic_regex_tokenizer import BasicRegexTokenizer

class Expression(BasicRegexTokenizer):
    def pattern(self):
        return r'\b((?<!over )expression|(?<!under )expression)\b'

    def token_type(self):
        return 'Expression'
