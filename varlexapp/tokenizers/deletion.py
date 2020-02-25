from .basic_regex_tokenizer import BasicRegexTokenizer

class Deletion(BasicRegexTokenizer):
    def pattern(self) -> str:
        return r'\b(del|deletion|(copy number)? ?loss)\b'

    def token_type(self) -> str:
        return 'Deletion'
