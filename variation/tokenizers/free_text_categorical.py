"""A module for free text categorical variation tokenization"""
from typing import Optional

from variation.schemas.token_response_schema import AmplificationToken
from variation.tokenizers.tokenizer import Tokenizer


class FreeTextCategorical(Tokenizer):
    """The Free Text Categorical tokenizer class"""

    def match(self, input_string: str) -> Optional[AmplificationToken]:
        """Return tokens that match the input string.
        Only supports amplification for now

        :param input_string: Input string
        :return: AmplificationToken token if a match is found
        """
        if input_string.lower() == "amplification":
            return AmplificationToken(token=input_string, input_string=input_string)

        return None
