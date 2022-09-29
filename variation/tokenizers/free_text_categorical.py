"""A module for free text categorical variation tokenization"""
from typing import Optional

from variation.schemas.token_response_schema import TokenMatchType, AmplificationToken,\
    Nomenclature
from variation.tokenizers import Tokenizer


class FreeTextCategorical(Tokenizer):
    """The Free Text Categorical tokenizer class"""

    def match(self, input_string: str) -> Optional[AmplificationToken]:
        """Return tokens that match the input string.
        Only supports amplification for now

        :param str input_string: Input string
        :return: AmplificationToken token if a match is found
        """
        if input_string.lower().strip() == "amplification":
            return AmplificationToken(
                token=input_string,
                input_string=input_string,
                match_type=TokenMatchType.UNSPECIFIED.value,
                nomenclature=Nomenclature.FREE_TEXT.value
            )
        return None
