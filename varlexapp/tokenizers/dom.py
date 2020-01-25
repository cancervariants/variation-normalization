from .basic_regex_tokenizer import BasicRegexTokenizer

class DOM(BasicRegexTokenizer):
    def pattern(self):
        return r'\b(\w+ )?domain\b'

    def token_type(self):
        return 'DOM'
