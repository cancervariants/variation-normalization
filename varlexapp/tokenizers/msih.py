from .basic_regex_tokenizer import BasicRegexTokenizer

class MSIH(BasicRegexTokenizer):
    def pattern(self):
        return r'\bMicrosatellite Instability-High\b'

    def token_type(self):
        return 'MSIH'
