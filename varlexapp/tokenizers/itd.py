from .basic_regex_tokenizer import BasicRegexTokenizer

class ITD(BasicRegexTokenizer):
    def pattern(self):
        return r'\b(Internal tandem duplication)|(ITD)\b'

    def token_type(self):
        return 'ITD'
