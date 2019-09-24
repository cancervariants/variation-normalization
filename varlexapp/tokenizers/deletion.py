from .basic_regex_tokenizer import BasicRegexTokenizer

class Deletion(BasicRegexTokenizer):
    def pattern(self):
        return r'\b(del|deletion|(copy number)? ?loss)\b'

    def token_type(self):
        return 'Deletion'
