from .basic_regex_tokenizer import BasicRegexTokenizer

class Methylation(BasicRegexTokenizer):
    def pattern(self):
        return r'(METHYLATION|HYPERMETHYLATION|HYPOMETHYLATION)'

    def token_type(self):
        return 'methylation'
