from .basic_regex_tokenizer import BasicRegexTokenizer

class LOCAL(BasicRegexTokenizer):
    def pattern(self):
        return r'\b(NUCLEAR EXPRESSION)|(CYTOPLASMIC EXPRESSION)|(CYTOPLASMIC MISLOCALIZATION)\b'

    def token_type(self):
        return 'LOCAL'
