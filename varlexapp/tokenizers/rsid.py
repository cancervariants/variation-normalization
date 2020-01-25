from .basic_regex_tokenizer import BasicRegexTokenizer

class RSID(BasicRegexTokenizer):
    def pattern(self):
        return r'\b(RS|rs|SNP|snp)\d*\s*\w*\b'

    def token_type(self):
        return 'RSID'
