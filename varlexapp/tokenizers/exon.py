from .basic_regex_tokenizer import BasicRegexTokenizer

class Exon(BasicRegexTokenizer):
    def pattern(self):
        return r'\bexon \d+\b'

    def token_type(self):
        return 'Exon'
