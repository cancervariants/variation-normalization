from .basic_regex_tokenizer import BasicRegexTokenizer

class WildType(BasicRegexTokenizer):
    def pattern(self):
        return r'\b(wild type)|wt\b'

    def token_type(self):
        return 'WildType'
