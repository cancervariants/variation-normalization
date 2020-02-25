from .basic_regex_tokenizer import BasicRegexTokenizer

class Amplification(BasicRegexTokenizer):
    def pattern(self) -> str:
        return r'\b(amp|amplification|(copy number)? ?gain)\b'

    def token_type(self) -> str:
        return 'Amplification'
