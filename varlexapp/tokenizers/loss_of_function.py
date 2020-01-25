from .basic_regex_tokenizer import BasicRegexTokenizer

class LossOfFunction(BasicRegexTokenizer):
    def pattern(self):
        return r'\b(LOSS[ -]OF[ -]FUNCTION)|(BIALLELIC INACTIVATION)\b'

    def token_type(self):
        return 'LossOfFunction'
