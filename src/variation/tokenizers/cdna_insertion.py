"""A module for Cdna Insertion Tokenization."""

from cool_seq_tool.schemas import AnnotationLayer

from variation.regex import CDNA_GENOMIC_INSERTION
from variation.schemas.token_response_schema import CdnaInsertionToken
from variation.tokenizers.tokenizer import Tokenizer


class CdnaInsertion(Tokenizer):
    """Class for tokenizing Insertion at the cdna reference sequence."""

    def match(self, input_string: str) -> CdnaInsertionToken | None:
        """Return a CdnaInsertionToken match if one exists.

        :param input_string: The input string to match
        :return: A CdnaInsertionToken if a match exists.
            Otherwise, None.
        """
        og_input_string = input_string
        _, input_string = self.strip_coord_prefix(
            input_string, match_coord_type=AnnotationLayer.CDNA
        )
        if not input_string:
            return None

        match = CDNA_GENOMIC_INSERTION.match(input_string)

        if match:
            match_dict = match.groupdict()
            pos0 = int(match_dict["pos0"])
            pos1 = int(match_dict["pos1"])
            inserted_sequence = match_dict["inserted_sequence"]

            return CdnaInsertionToken(
                input_string=og_input_string,
                token=f"{pos0}_{pos1}{inserted_sequence}",
                pos0=pos0,
                pos1=pos1,
                inserted_sequence=inserted_sequence,
            )

        return None
