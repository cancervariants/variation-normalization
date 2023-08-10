"""A module for Genomic Duplication Tokenization."""
from typing import Optional

from cool_seq_tool.schemas import AnnotationLayer

from variation.regex import (
    GENOMIC_DUPLICATION,
    GENOMIC_DUPLICATION_AMBIGUOUS_1,
    GENOMIC_DUPLICATION_AMBIGUOUS_2,
    GENOMIC_DUPLICATION_AMBIGUOUS_3,
)
from variation.schemas.app_schemas import AmbiguousRegexType
from variation.schemas.token_response_schema import (
    GenomicDuplicationAmbiguousToken,
    GenomicDuplicationToken,
)
from variation.tokenizers.tokenizer import Tokenizer


class GenomicDuplication(Tokenizer):
    """Class for tokenizing duplications on the genomic coordinate."""

    def match(self, input_string: str) -> Optional[GenomicDuplicationToken]:
        """Return a GenomicDelInsToken match if one exists.

        :param input_string: The input string to match
        :return: A GenomicDelInsToken if a match exists.
            Otherwise, None.
        """
        og_input_string = input_string
        _, input_string = self.strip_coord_prefix(
            input_string, match_coord_type=AnnotationLayer.GENOMIC
        )
        if not input_string:
            return None

        # First try matching on simple genomic duplications
        match = GENOMIC_DUPLICATION.match(input_string)

        if match:
            match_dict = match.groupdict()

            return GenomicDuplicationToken(
                input_string=og_input_string,
                token=input_string,
                pos0=int(match_dict["pos0"]),
                pos1=int(match_dict["pos1"]) if match_dict["pos1"] else None,
            )
        else:
            # Going to try ambiguous genomic duplications
            match = GENOMIC_DUPLICATION_AMBIGUOUS_1.match(input_string)
            if match:
                match_dict = match.groupdict()
                pos0 = match_dict["pos0"]
                pos1 = match_dict["pos1"]
                pos2 = match_dict["pos2"]
                pos3 = match_dict["pos3"]

                # (?_?)_(#_#), (#_#)_(?, ?), (?_?)_(?_?) are not supported
                if not any(
                    ((pos0 == "?" and pos1 == "?"), (pos2 == "?" and pos3 == "?"))
                ):
                    return GenomicDuplicationAmbiguousToken(
                        input_string=og_input_string,
                        token=input_string,
                        pos0=int(pos0) if pos0 != "?" else pos0,
                        pos1=int(pos1) if pos1 != "?" else pos1,
                        pos2=int(pos2) if pos2 != "?" else pos2,
                        pos3=int(pos3) if pos3 != "?" else pos3,
                        ambiguous_regex_type=AmbiguousRegexType.REGEX_1,
                    )

            else:
                for pattern_re, regex_type in [
                    (GENOMIC_DUPLICATION_AMBIGUOUS_2, AmbiguousRegexType.REGEX_2),
                    (GENOMIC_DUPLICATION_AMBIGUOUS_3, AmbiguousRegexType.REGEX_3),
                ]:
                    match = pattern_re.match(input_string)

                    if match:
                        matched_pos = dict()
                        match_dict = match.groupdict()
                        for k in match_dict:
                            v = match_dict[k]
                            if v:
                                v = int(v) if v != "?" else v

                            matched_pos[k] = v

                        return GenomicDuplicationAmbiguousToken(
                            input_string=og_input_string,
                            token=input_string,
                            pos0=matched_pos["pos0"],
                            pos1=matched_pos.get("pos1"),
                            pos2=matched_pos["pos2"],
                            pos3=matched_pos.get("pos3"),
                            ambiguous_regex_type=regex_type,
                        )
